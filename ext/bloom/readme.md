## Description
**C++ Bloom Filter Library**, has the following capabilities:

+ Optimal parameter selection based on expected false positive rate.
+ Union, intersection and difference operations between bloom filters.
+ Compression of in-use table (increase of false positive probability vs space)
+ Portable and efficient source code implementation.
+ Single header implementation, no building required. No external dependencies


## Compatible Compilers
+ GNU Compiler Collection (4.1+)
+ Intel® C++ Compiler (9.x+)
+ Clang/LLVM (1.1+)
+ PGI C++ (10.x+)
+ Microsoft Visual Studio C++ Compiler (8.1+)

For more information please visit: http://www.partow.net/programming/bloomfilter/index.html

---

## Simple Bloom Filter Example
**Example's objectives:**
+ Instantiate and configure a Bloom filter
+ Add some strings and integers to the Bloom filter
+ Query the Bloom filter for membership of the previously added strings and integers
+ Query the Bloom filter for membership of integers that were **NOT** previously added *(potential false positives)*

```javascript
#include <iostream>
#include <string>

#include "bloom_filter.hpp"

int main()
{

   bloom_parameters parameters;

   // How many elements roughly do we expect to insert?
   parameters.projected_element_count = 1000;

   // Maximum tolerable false positive probability? (0,1)
   parameters.false_positive_probability = 0.0001; // 1 in 10000

   // Simple randomizer (optional)
   parameters.random_seed = 0xA5A5A5A5;

   if (!parameters)
   {
      std::cout << "Error - Invalid set of bloom filter parameters!" << std::endl;
      return 1;
   }

   parameters.compute_optimal_parameters();

   //Instantiate Bloom Filter
   bloom_filter filter(parameters);

   std::string str_list[] = { "AbC", "iJk", "XYZ" };

   // Insert into Bloom Filter
   {
      // Insert some strings
      for (std::size_t i = 0; i < (sizeof(str_list) / sizeof(std::string)); ++i)
      {
         filter.insert(str_list[i]);
      }

      // Insert some numbers
      for (std::size_t i = 0; i < 100; ++i)
      {
         filter.insert(i);
      }
   }


   // Query Bloom Filter
   {
      // Query the existence of strings
      for (std::size_t i = 0; i < (sizeof(str_list) / sizeof(std::string)); ++i)
      {
         if (filter.contains(str_list[i]))
         {
            std::cout << "BF contains: " << str_list[i] << std::endl;
         }
      }

      // Query the existence of numbers
      for (std::size_t i = 0; i < 100; ++i)
      {
         if (filter.contains(i))
         {
            std::cout << "BF contains: " << i << std::endl;
         }
      }

      // Query the existence of invalid numbers
      for (int i = -1; i > -100; --i)
      {
         if (filter.contains(i))
         {
            std::cout << "BF falsely contains: " << i << std::endl;
         }
      }
   }

   return 0;
}
```

