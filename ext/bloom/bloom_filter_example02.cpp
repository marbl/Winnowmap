/*
 *********************************************************************
 *                                                                   *
 *                           Open Bloom Filter                       *
 *                                                                   *
 * Description: Demonstration of a Bloom Filter                      *
 * Author: Arash Partow - 2000                                       *
 * URL: http://www.partow.net                                        *
 * URL: http://www.partow.net/programming/hashfunctions/index.html   *
 *                                                                   *
 * Copyright notice:                                                 *
 * Free use of the Open Bloom Filter Library is permitted under the  *
 * guidelines and in accordance with the MIT License.                *
 * http://www.opensource.org/licenses/MIT                            *
 *                                                                   *
 *********************************************************************
*/


/*
   Description: This example will demonstrate how to instantiate a Bloom filter,
                insert strings and then query the inserted strings and a set of
                outlier strings for membership status within the Bloom filter.
                Furthermore this process will be carried out upon 1000 unique
                instances of Bloom filter. The objective is to empirically
                determine which "random" seed that when used to construct a
                Bloom filter will provide the smallest observed false positive
                probability for the given sets of data. The optimal seed will
                be the one associated with the round that has the smallest
                difference percentage of false positive probability against
                the user specified false positive probability.
*/


#include <iostream>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <deque>
#include <set>
#include <string>

#include "bloom_filter.hpp"

bool load_word_list(int argc, char* argv[], std::vector<std::string>& word_list);

template <class T,
          class Allocator,
          template <class,class> class Container>
bool read_file(const std::string& file_name, Container<T, Allocator>& c);

std::string uppercase(std::string str);
std::string reverse(std::string str);

void generate_outliers(const std::vector<std::string>& word_list, std::deque<std::string>& outliers);
void purify_outliers(const std::vector<std::string>& word_list,std::deque<std::string>& outliers);

int main(int argc, char* argv[])
{
   std::vector<std::string> word_list;
   std::deque<std::string> outliers;

   if (!load_word_list(argc,argv,word_list))
   {
      return 1;
   }

   generate_outliers(word_list,outliers);

   unsigned int random_seed           = 0;
   std::size_t word_list_storage_size = 0;

   for (unsigned int i = 0; i < word_list.size(); ++i)
   {
      word_list_storage_size += word_list[i].size();
   }

   std::size_t total_number_of_queries = 0;

   const double desired_probability_of_false_positive = 1.0 / word_list.size();

   printf("Round\t   Queries\t   FPQ\t   IPFP\t           PFP\t            DPFP\t    TvD\n");

   unsigned int max_false_positive_count = 0;
   unsigned int min_false_positive_count = std::numeric_limits<unsigned int>::max();
   unsigned int total_false_positive     = 0;
   unsigned int total_zero_fp            = 0;
   unsigned long long int bloom_filter_size = 0;


   static const unsigned int rounds = 1000;

   while (random_seed < rounds)
   {
      bloom_parameters parameters;

      parameters.projected_element_count    = word_list.size();
      parameters.false_positive_probability = desired_probability_of_false_positive;
      parameters.random_seed                = ++random_seed;

      if (!parameters)
      {
         std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
         return 1;
      }

      parameters.compute_optimal_parameters();

      bloom_filter filter(parameters);

      filter.insert(word_list.begin(),word_list.end());

      std::vector<std::string>::iterator it = filter.contains_all(word_list.begin(),word_list.end());

      if (word_list.end() != it)
      {
         std::cout << "ERROR: key not found in bloom filter! =>" << (*it) << std::endl;
         return 1;
      }

      unsigned int current_total_false_positive = 0;

      for (std::deque<std::string>::iterator itr = outliers.begin(); itr != outliers.end(); ++itr)
      {
         if (filter.contains(*itr))
         {
            ++current_total_false_positive;
         }
      }

      total_number_of_queries += (outliers.size() + word_list.size());

      // Overall false positive probability
      double pfp = current_total_false_positive / (1.0 * outliers.size());

      printf("%6llu\t%10llu\t%6d\t%8.7f\t%8.7f\t%9.3f%%\t%8.6f\n",
             static_cast<unsigned long long>(random_seed),
             static_cast<unsigned long long>(total_number_of_queries),
             current_total_false_positive,
             desired_probability_of_false_positive,
             pfp,
             (100.0 * pfp) / desired_probability_of_false_positive,
             (100.0 * filter.size()) / (bits_per_char * word_list_storage_size));

      if (current_total_false_positive < min_false_positive_count)
         min_false_positive_count = current_total_false_positive;
      else if (current_total_false_positive > max_false_positive_count)
         max_false_positive_count = current_total_false_positive;

     total_false_positive += current_total_false_positive;

     if (0 == current_total_false_positive)
        ++total_zero_fp;

     bloom_filter_size = filter.size();
   }

   double average_fpc = (1.0 * total_false_positive) / rounds;
   double average_fpp = average_fpc / (outliers.size() + word_list.size());

   printf("Bloom Filter Statistics\n"
          "MinFPC: %d\tMaxFPC: %d\tAverageFPC: %8.5f\tAverageFPP: %9.8f Zero-FPC:%d\n"
          "Filter Size: %lluKB\tData Size: %dKB\n",
          min_false_positive_count,
          max_false_positive_count,
          average_fpc,
          average_fpp,
          total_zero_fp,
          bloom_filter_size / (8 * 1024),
          static_cast<unsigned int>(word_list_storage_size / 1024));

   /*
      Terminology
      MinFPC : Minimum (smallest) False Positive Count
      MaxFPC : Maximum (largest) False Positive Count
      AverageFPC : Average False Positive Count
      AverageFPP : Average False Positive Probability

      FPQ  : False Positive Queries
      IPFP : Indicative (desired) False Positive Probability
      PFP  : Probability of a False Positive (based on the FPQ)
      DPFP : Difference as a percentage between IPFP and PFP
      TvD  : percentage of the filter size versus the raw data size
   */

   return 0;
}

bool load_word_list(int argc, char* argv[], std::vector<std::string>& word_list)
{
   // Note: The word-lists can be obtained from:
   // https://github.com/ArashPartow/bloom
   static const std::string wl_list[] =
                     { "word-list.txt",
                       "word-list-large.txt",
                       "word-list-extra-large.txt",
                       "random-list.txt"
                     };

   std::size_t index = 0;

   if (2 == argc)
   {
      index = ::atoi(argv[1]);

      const std::size_t wl_list_size = sizeof(wl_list) / sizeof(std::string);

      if (index >= wl_list_size)
      {
         std::cout << "Invalid world list index: " << index << std::endl;
         return false;
      }
   }

   std::cout << "Loading list " << wl_list[index] << ".....";

   if (!read_file(wl_list[index],word_list))
   {
      return false;
   }

   if (word_list.empty())
   {
      std::cout << "No word list - Either none requested, or desired word list could not be loaded." << std::endl;
      return false;
   }
   else
      std::cout << " Complete." << std::endl;

   return true;
}

template <class T,
          class Allocator,
          template <class,class> class Container>
bool read_file(const std::string& file_name, Container<T, Allocator>& c)
{
   std::ifstream stream(file_name.c_str());

   if (!stream)
   {
      std::cout << "Error: Failed to open file '" << file_name << "'" << std::endl;
      return false;
   }

   std::string buffer;

   while (std::getline(stream,buffer))
   {
      c.push_back(buffer);
      c.push_back(uppercase(buffer));
   }

   return true;
}

std::string uppercase(std::string str)
{
   for (std::size_t i = 0; i < str.size(); ++i)
   {
      str[i] = static_cast<char>(toupper(str[i]));
   }

   return str;
}

std::string reverse(std::string str)
{
   // Not the most efficient way of doing this.
   std::reverse(str.begin(),str.end());
   return str;
}

void generate_outliers(const std::vector<std::string>& word_list, std::deque<std::string>& outliers)
{
   std::cout << "Generating outliers..... ";

   for (std::vector<std::string>::const_iterator it = word_list.begin(); it != word_list.end(); ++it)
   {
      if ((*it) != reverse((*it)))
      {
         outliers.push_back((*it) + reverse((*it)));
         outliers.push_back((*it) + (*it));
         outliers.push_back(reverse((*it)) + (*it) + reverse((*it)));
      }

      std::string ns = *it;

      for (unsigned int i = 0; i < ns.size(); ++i)
      {
         if (1 == (i & 0x00)) ns[i] = ~ns[i];
      }

      outliers.push_back(ns);
   }

   static const std::string rand_str[] =
                  {
                     "oD5l", "pccW", "5yHt", "ndaN", "OaJh", "tWPc", "Cr9C", "a9zE",
                     "H1wL", "yo1V", "16D7", "f2WR", "0MVQ", "PkKn", "PlVa", "MvzL",
                     "9Csl", "JQTv", "IveD", "FDVS", "Q7HE", "QgcF", "Q9Vo", "V8zJ",
                     "EJWT", "GuLC", "rM3d", "PJF4", "HXPW", "qKx3", "ztRP", "t4KP",
                     "m1zV", "fn12", "B1QP", "Jr4I", "Mf8M", "4jBd", "anGR", "Pipt",
                     "QHon", "GNlc", "UeXM", "mVM5", "ABI8", "RhB3", "5h2s", "hOYo",
                     "gaId", "DX40", "THMu", "EwlP", "n9Mz", "oC1S", "BfMl", "uCZ1",
                     "G2bA", "MOH9", "zZ0O", "PKDO", "3nRU", "Z6ie", "4cso", "LnQO",
                     "MJTtT","td3rC","A5JNR","1yL5B","rQnJk","jNKYF","CD0XD","pFLSG",
                     "fxO1a","CAjBE","ORk4e","0LERI","R7d0x","Qqd7v","6Kih5","9tTCB",
                     "yCg9U","D2Tv7","XpNHn","6zeFQ","BT2cs","WGhKW","zTv6B","TTPFk",
                     "XjNVX","pg9yW","4pKiZ","mQUhL","xrXzR","kVRm5","NSyC4","olXm9",
                     "UWkYy","8Ys6r","yd4Fl","5L4mB","nP3nH","f0DFb","glnQa","DlXQa",
                     "cQdH6","eBmIN","fDj6F","ezLow","C15vu","I2Z2j","BQgzg","eVBid",
                     "hn5TO","WZyQN","xXgsE","sL6nK","8DKD8","jcrbp","AcRak","h8N5o",
                     "LViwC","ThEKf","O7fd5","oN0Id","OM1m0","4OLiR","VIa8N","bJZFG",
                     "9j3rL","SzW0N","7m7pY","mY9bg","k1p3e","3OFm1","r45se","VYwz3",
                     "pDjXt","ZcqcJ","npPHx","hA3bw","w7lSO","jEmZL","1x3AZ","FN47G",
                     "kThNf","aC4fq","rzDwi","CYRNG","gCeuG","wCVqO","d1R60","bEauW",
                     "KeUwW","lIKhO","RfPv3","dK5wE","1X7qu","tRwEn","1c03P","GwHCl",
                     "CsJaO","zl4j1","e0aEc","Uskgi","rgTGR","jyR4g","Tt6l4","lRoaw",
                     "94ult","qZwBX","eYW8S","Qf6UH","AbV56","N1hJq","JIaVe","8LHEx",
                     "DeNbS","30I0a","hm6qw","3jcaO","4WkuA","mQ219","Gb81C","yx4HM",
                     "Chlqfi9S1y", "BMwUgVFu2X", "ZmpEGOVrVe", "13ggJxrPkC", "fcJJpyMGjm", "9T00Dv4ZAb",
                     "p3YRcP7M2o", "sR0qNUXCHv", "gCxWZbJ6rb", "R4YtzRXXUl", "vwyYz5j6pY", "XPWUvLXhJ7",
                     "7PwfnVVb7U", "1f34Q6hOYz", "1EM2abZY61", "0a6Ivi4S0a", "Teq2LrQs2T", "dWXLCgWHc8",
                     "LawMv7ujn4", "N8VFgbZQx5", "tfvHHxoDgi", "ImwYgXA2tf", "KkIES9NqZO", "ajcz0qjjda",
                     "6Vz28vlGs9", "VMCc5W8cCt", "BiQB8BRJ98", "43CpOJSMpA", "jfBJdqwXcU", "ecHR9EO2ib",
                     "LH7CcXyCZ7", "JntqGSgSpa", "0MbTMpZPFW", "5FJSdiCXzR", "5gda2AhA2x", "lrDFc1lnXk",
                     "zrEwECHvjs", "B0JldDxFa1", "6DYal4QxKa", "Hsqx6kP2S4", "zZwnALSuFh", "Shh4ISZcKW",
                     "P9VDaNSk7Z", "mEI2PLSCO6", "WyTyrQORtu", "IvJyMMRgh3", "Q6pgJq8Nkv", "dhOgR3tDAD",
                     "Y9h6bVgbxO", "wA15tiOPTm", "8TaIKf1zCO", "z75dzabHBs", "AS6OPnwoJI", "2DSZka9Auj",
                     "QLzUjV2CWs", "KZSN2SVhia", "7ttYKWF2ue", "1Zxfu7B2ST", "RnkpmwjsCi", "YpcSIzaqx5",
                     "RDEwFD9gmX", "Nlx3V4Cjw4", "9ZdvITOj8M", "httUPWMNXO", "Ypv9PjxGwa", "LlwyNolNnH",
                     "6xpJOht47a", "tbmz4WIdcG", "OwzuVDlb7D", "PBQKJxo8DQ", "uVnMQn7hK6", "rlnZINuDUa",
                     "2feyyYukPa", "teOlpKuDBn", "LxBSWh0dL1", "Onyb7r4Jp0", "bZxXE6xOXg", "d9NSvNTunQ",
                     "ONerLBic32", "8mar4rKmFk", "5cCN9uwaCg", "ElVrYOHHMv", "YF6Og8DX40", "OgiCwpCQ5a",
                     "K6nSRZVxdR", "gqyXXXoVFW", "ulyRYizcBP", "khUx31K5UR", "qZFRzVthju", "pQBh0vnB20",
                     "dk8NIN7ajy", "XP7ed1OjZx", "IRYNwA5iFR", "hiSEBhTukC", "Ns4jJ3jzGo", "dYoCSxjIvM",
                     "HzGLbl5i1g", "baizENd4ko", "6rCqGBO8t1", "QWGfC8UaA7", "JFhRfxQe4K", "8R4W6IWANz",
                     "2TnWf1w7JH", "0z69e0wcoG", "8SN1mRHCY7", "oFGCYHHwGX", "G8xqnBgxjO", "6B3SAOayHt",
                     "XRW3ZSG1gw", "WcIjTxMxOM", "wNqCAIaTb2", "gO4em4HW8H", "TgGFSMEtbG", "WiwmbEw3QA",
                     "D2xshYUgpu", "xRUZCQVzBs", "nCnUmMgIjE", "p4Ewt1yCJr", "MeOjDcaMY5", "1XelMeXiiI"
                  };

   static const std::size_t rand_str_size = sizeof(rand_str) / sizeof(std::string);

   for (unsigned int i = 0; i < rand_str_size; ++i)
   {
      std::string s0 = rand_str[i];
      std::string s1 = rand_str[(i + 1) % rand_str_size];
      std::string s2 = rand_str[(i + 2) % rand_str_size];
      std::string s3 = rand_str[(i + 3) % rand_str_size];
      std::string s4 = rand_str[(i + 4) % rand_str_size];
      std::string s5 = rand_str[(i + 5) % rand_str_size];
      std::string s6 = rand_str[(i + 6) % rand_str_size];

      outliers.push_back(s0);
      outliers.push_back(s0 + s1);
      outliers.push_back(s0 + s2 + s4);
      outliers.push_back(s0 + s1 + s3);
      outliers.push_back(s0 + s1 + s2 + s3 + s4 + s5);
      outliers.push_back(s0 + s1 + s2 + s3 + s4 + s5 + s6);

      outliers.push_back(reverse(s0));
      outliers.push_back(reverse(s0 + s1));
      outliers.push_back(reverse(s0 + s2 + s4));
      outliers.push_back(reverse(s0 + s1 + s3));
      outliers.push_back(reverse(s0 + s1 + s2 + s3 + s4 + s5));
      outliers.push_back(reverse(s0 + s1 + s2 + s3 + s4 + s5 + s6));
   }

   std::sort(outliers.begin(),outliers.end());

   purify_outliers(word_list,outliers);

   std::cout << "Complete." << std::endl;
}

void purify_outliers(const std::vector<std::string>& word_list, std::deque<std::string>& outliers)
{
   std::set<std::string> set1;
   std::set<std::string> set2;

   std::copy(word_list.begin(), word_list.end(),std::inserter(set1,set1.begin()));
   std::copy(outliers.begin(), outliers.end(), std::inserter(set2,set2.begin()));

   std::deque<std::string> intersect_list;

   std::set_intersection(set1.begin(),set1.end(),
                         set2.begin(),set2.end(),
                         std::back_inserter(intersect_list));

   std::sort(intersect_list.begin(),intersect_list.end());

   if (!intersect_list.empty())
   {
      std::deque<std::string> new_outliers;

      for (std::deque<std::string>::iterator it = outliers.begin(); it != outliers.end(); ++it)
      {
         if (!std::binary_search(intersect_list.begin(),intersect_list.end(),*it))
         {
            new_outliers.push_back(*it);
         }
      }

      outliers.swap(new_outliers);
   }
}
