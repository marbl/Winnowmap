#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of meryl-utility, a collection of miscellaneous code
 #  used by Meryl, Canu and others.
 #
 #  This software is based on:
 #    'Canu' v2.0              (https://github.com/marbl/canu)
 #  which is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

use strict;
use File::Compare;
use Cwd qw(getcwd);

my $cwd = getcwd();

#
#  Update version.H with module and version info.
#
#  Expects four parameters on the command line:
#   - module-name:        name of the module, e.g., 'verkko' or 'canu'.
#   - release-type:       type of release this is:
#                           'snapshot' - query git to discover the version number
#                           'release'  - use the next arg as the version number
#   - version-number:     the version string to report, e.g., v2.0.1
#                            MUST be of form 'vMAJOR.MINOR[.PATCH]'
#   - path-to-version.H:  path to file 'version.H'; 
#
#  The script will create output file 'path-to-version.H' and write
#  configuration information suitable for inclusion in both C source code and
#  Makefiles.
#
#  IMPORTANT: The script must output ONLY the full path to the output file.
#  Makefile uses this.  Any other output, EVEN to STDERR, will corrupt this.
#
#  If any errors occur, just print a bogus filename, and gmake will blow up.
#

my $modName  =  $ARGV[0];   #  Name of module we're making versions for.
my $MODNAME  =  $ARGV[0] =~ tr/a-z-/A-Z_/r;    #  And same in uppercase.
my $typein   =  $ARGV[1];   #  'snapshot' or 'release'.
my $versin   =  $ARGV[2];   #  'vM.N[.P]'.
my $major    =  "M";
my $minor    =  "N";
my $verFile  =  $ARGV[3];   #  Path to output version.H.

my $branch   = "master";
my $commits  = undef;       #  Local and output variables for creating
my $hash     = undef;       #  a version string using git.
my $revCount = 0;
my $dirty    = undef;
my @submodules;

if (($modName eq "") || ($modName eq undef) || ($verFile eq undef)) {
    die "usage: $0 module-name release-type version-number version-file.H\n";
}

#  Possibilities:
#    'release'  tag and no .git - an actual release!
#    'release'  tag and    .git - a maintenance build.
#    'snapshot' tag and no .git - an ERROR!
#    'snapshot' tag and    .git - a development build.

my $isgit = (system("git rev-parse --is-inside-work-tree > /dev/null 2>&1") == 0);


if    (($typein eq "release") && ($isgit == 0)) {   #  release tag, no git repo
  if ($versin =~ m!^v(\d+)\.(\d+\.?\d*)$!) {
    $branch = "release";
    $major  = $1;
    $minor  = $2;
  }
  else {
    die "Invalid version string '$typein' '$versin'.\n";
  }
}

elsif (($typein eq "release") && ($isgit != 0)) {   #  release tag and a git repo.
  countRevisions();
  findCommitVersion("release");
  findLocalChanges();
  findBranch();
  findSubmodules();
}

elsif (($typein ne "release") && ($isgit == 0)) {   #  development, but outside git?
  die "Development copy outside git?\n";
}

elsif (($typein ne "release") && ($isgit != 0)) {   #  developtment and a git repo.
  countRevisions();
  findCommitVersion("branch");
  findLocalChanges();
  findBranch("master");
  findSubmodules();
}

else {
  die "Logic failed.\n";
}

updateVersionFile();
exit(0);




#
#  FUNCTIONS BELOW
#


#
#  Count the number of changes since the last release.
#
#  We used to also remember the first hash returned (as $hash2) but it
#  isn't used anymore.
#
sub countRevisions () {
  $revCount = 0;

  open(F, "git rev-list HEAD |") or die "Failed to run 'git rev-list'.\n";
  while (<F>) {
    $revCount++;
  }
  close(F);
}

#
#  Find the commit and version we're at.
#
sub findCommitVersion () {
  open(F, "git describe --tags --long --always --abbrev=40 |") or die "Failed to run 'git describe'.\n";
  while (<F>) {
    if (m/^v(\d+)\.(\d+).*-(\d+)-g(.{40})$/)  {
      ($major, $minor, $commits, $hash) = ($1, $2, $3, $4);
    }
  }
  close(F);
}

#
#  Decide if we've got locally committed changes.
#
sub findLocalChanges () {
  my ($da, $dc);

  open(F, "git status |");
  while (<F>) {
    if (m/is\s+ahead\s/)                 { $da = "ahead of github";  }
    if (m/not\s+staged\s+for\s+commit/)  { $dc = "w/changes";        }
  }
  close(F);
  $da = undef;
  $dc = undef;

  if    (defined($da) && defined($dc)) { $dirty = "ahead of github w/changes"; }
  elsif (defined($da))                 { $dirty = "ahead of github";           }
  elsif (defined($dc))                 { $dirty = "w/changes";                 }
  else                                 { $dirty = "sync'd with github";        }
}

#
#  Figure out which branch we're on.
#
sub findBranch ($) {
  open(F, "git branch --contains |");   #  Show branches that contain the tip commit.
  while (<F>) {
    next  if (m/detached\sat/);         #  Skip 'detached' states.

    s/^\*{0,1}\s+//;                    #  Remove any '*' annotation
    s/\s+$//;                           #  and all spaces.

    $branch = $_;                       #  Pick the first branch
    last;                               #  mentioned.
  }
  close(F);

  if ($branch =~ m/v(\d+)\.(\d+)/) {    #  If a release branch,
    $major = $1;                        #  save major and minor
    $minor = $2;                        #  version numbers.
  }
}

#
#  Get information on any submodules present.
#
sub findSubmodules () {
  open(F, "git submodule status |");
  while (<F>) {
    if (m/^(.*)\s+(.*)\s+\((.*)\)$/) {
      push @submodules, "$2 $3 $1";
    }
  }
  close(F);
}



#
#  Create a version file, but don't change the timestamp if there
#  are no changes (otherwise, we rebuild everything on every build!).
#
sub updateVersionFile () {
  my $cv;
  my $mv;

  if (defined($commits))  { $cv = "$branch +$commits changes (r$revCount $hash)"; }
  else                    { $cv = "$major.$minor"; }

  if (defined($commits))  { $mv = "$branch +$commits changes (r$revCount $hash) ($dirty)"; }
  else                    { $mv = "$branch v$major.$minor"; }

  $cv =~ s/\+1\schanges/+1 change/;  #  Easier to fix up than
  $mv =~ s/\+1\schanges/+1 change/;  #  to do correctly above.

  open(F, "> $verFile.new") or die "Failed to open '$verFile.new' for writing: $!\n";
  print F "#if 0\n";
  print F "#/*\n";
  print F "#  Automagically generated by version_update.pl!  Do not commit!\n";
  print F "#*/\n";
  print F "#endif\n";
  print F "\n";
  print F "#if 0\n";
  print F "#/*\n";
  print F "#  These lines communicate the version information to C++.\n";
  print F "#*/\n";
  print F "#endif\n";
  print F "#define ${MODNAME}_VERSION_BRANCH    \"$branch\"\n";
  print F "#define ${MODNAME}_VERSION_MAJOR     \"$major\"\n";
  print F "#define ${MODNAME}_VERSION_MINOR     \"$minor\"\n";
  print F "#define ${MODNAME}_VERSION_COMMITS   \"$commits\"\n";
  print F "#define ${MODNAME}_VERSION_REVISION  \"$revCount\"\n";
  print F "#define ${MODNAME}_VERSION_HASH      \"$hash\"\n";          #  SAM/BAM support REQUIRES
  print F "\n";                                                        #  that there is no newline
  print F "#define ${MODNAME}_VERSION           \"$modName $cv\"\n";   #  on ${MODNAME}_VERSION.
  print F "\n"                                                   if ($MODNAME ne "MERYL_UTILITY");
  print F "#undef  MERYL_UTILITY_VERSION\n"                      if ($MODNAME ne "MERYL_UTILITY");
  print F "#define MERYL_UTILITY_VERSION ${MODNAME}_VERSION\n"   if ($MODNAME ne "MERYL_UTILITY");
  print F "\n";
  print F "#if 0\n";
  print F "#/*\n";
  print F "#  These lines communicate version information to make.\n";
  print F "define BUILDING\n";
  print F "Building $mv\n";
  print F "  with $_\n"     foreach (@submodules);
  print F "endef\n";
  print F "#*/\n";
  print F "#endif\n";
  close(F);

  if (compare("$verFile", "$verFile.new") == 0) {    #  If they're the same, don't replace
    unlink "$verFile.new";                           #  the original.  This maintains
  } else {                                           #  the timestamp, which makes make
    unlink "$verFile";                               #  much happier.
    rename "$verFile.new", "$verFile";
  }

  print "$verFile\n";
}
