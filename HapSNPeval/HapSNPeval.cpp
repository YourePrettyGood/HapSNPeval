/*********************************************************************************
 * HapSNPeval.cpp                                                                *
 * Written by Patrick Reilly                                                     *
 * Version 1.0 written 2015/09/15                                                *
 * Description:                                                                  *
 *  This program evaluates the accuracy of reconstructed haplotypes as compared  *
 *  to the original (i.e. source for simulated reads) haplotypes.  The metrics   *
 *  of accuracy include switch count (number of times the reconstructed haplotype*
 *  switches identities at heterozygous SNPs), false SNP count, etc.             *
 *                                                                               *
 * Syntax: HapSNPeval -p true_haplotype_prefix input_alignment.fa                *
 *  input_alignment.fa:   Path to the multiple sequence alignment of the two true*
 *                        haplotypes with the two test haplotypes, in alignment  *
 *                        FASTA format                                           *
 *  true_haplotype_prefix:Prefix of the header string for each true haplotype    *
 *                                                                               *
 * Design:                                                                       *
 *  The general idea is to identify and read in the true and test haplotype      *
 *  sequences, then iterate along the two true haplotypes, identifying true      *
 *  heterozygous SNPs and assessing the identity of each test haplotype relative *
 *  to the previous identity of that test haplotype (to count switching errors), *
 *  as well as assessing deviations from homozygosity in the test haplotypes at  *
 *  homozygous sites in the true haplotypes, etc.                                *
 *  This program should allow high throughput assessment of haplotypes           *
 *  reconstructed from simulated read data using various haplotype assembly      *
 *  software/pipelines, primarily with regards to SNP phasing.                   *
 *********************************************************************************/
 
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
 
using namespace std;
 
int main(int argc, char *argv[]) {
   //Argument parsing variables:
   int helpflag = 0;
   int position_output_flag = 0;
   int optvalue;
   int optindex = 0;
   struct option long_options[] = 
      {
         {"help", no_argument, &helpflag, 1},
         {"position_output", no_argument, &position_output_flag, 1},
         {"true_prefix", required_argument, 0, 'p'},
         {0,0,0,0}
      };
   string true_prefix, input_alignment_file;
   //Core algorithm variables:
   string true_one_header = "", true_two_header = "", test_one_header = "", test_two_header = "", line_buffer;
   string true_one, true_two, test_one, test_two;
   unsigned short int test_one_id = 0, test_two_id = 0;
   unsigned long int test_one_switches = 0, test_two_switches = 0,
                     test_one_false_snps = 0, test_two_false_snps = 0,
                     test_one_false_indels = 0, test_two_false_indels = 0,
                     test_one_bad_calls = 0, test_two_bad_calls = 0;
   ifstream input_alignment;
   unsigned short int record_num = 0;
   
   //Parse input arguments with getopt_long:
   while ((optvalue = getopt_long(argc, argv, "hop:", long_options, &optindex)) != -1) {
      switch (optvalue) {
         case 0:
            //Flag was set, so skip
            break;
         case 'h':
            helpflag = 1;
            break;
         case 'o':
            position_output_flag = 1;
            break;
         case 'p':
            //Set the true haplotype prefix
            if (optarg == 0) {
               cerr << "Missing true haplotype prefix argument." << endl;
               helpflag = 3;
               break;
            }
            true_prefix = optarg;
            break;
         default:
            //Invalid argument
            cerr << "Invalid argument " << optvalue << " supplied." << endl;
            helpflag = 4;
            break;
      }
   }
   if (optind < argc) { //Read in the non-option argument, ignore any others
      input_alignment_file = argv[optind];
      input_alignment.open(input_alignment_file, ios_base::in);
      if (!input_alignment) {
         cerr << "Unable to open input alignment file." << endl;
         helpflag = 5;
      }
   } else { //Missing input alignment file path
      cerr << "Missing input alignment file path." << endl;
      helpflag = 6;
   }
   if (helpflag) { //If input errors or the help flag were detected, output usage and exit
      cout << "Usage: " << argv[0] << " -p true_haplotype_prefix input_alignment.fa" << endl;
      cout << " p\t\t\tPrefix of the header string for each true haplotype" << endl;
      cout << " input_alignment.fa\tPath to the MSA in FASTA format" << endl;
      return helpflag;
   }
   
   //WARNING: If too long of a haplotype is input using unwrapped FASTA, memory allocation issues may occur.
   //Could resolve this by using buffered binary reads, but for version 1.0 we will ignore it.
   //Read in the alignment records:
   while (input_alignment.good()) {
      getline(input_alignment, line_buffer);
      if (line_buffer[0] == '>') { //Header line
         if (line_buffer.find(true_prefix, 1) != string::npos) { //True haplotype record
            if (true_one_header == "") { //First true haplotype record
               true_one_header = line_buffer.substr(1);
               record_num = 1;
            } else { //Second true haplotype record
               true_two_header = line_buffer.substr(1);
               record_num = 2;
            }
         } else { //Test haplotype record
            if (test_one_header == "") { //First test haplotype record
               test_one_header = line_buffer.substr(1);
               record_num = 3;
            } else { //Second test haplotype record
               test_two_header = line_buffer.substr(1);
               record_num = 4;
            }
         }
      } else { //FASTA line
         //Since newlines are discarded, we can simply append each buffered line to the appropriate record
         switch (record_num) {
            case 1:
               true_one.append(line_buffer);
               break;
            case 2:
               true_two.append(line_buffer);
               break;
            case 3:
               test_one.append(line_buffer);
               break;
            case 4:
               test_two.append(line_buffer);
               break;
            default:
               break;
         }
         //Note: Extra newlines at the end of the FASTA file are handled (discarded during getline, so append operates on "").
      }
   }
   if (!input_alignment.eof()) { //Loop was not exited on EOF, so an error occurred
      cerr << "An error occurred while reading the input alignment file." << endl;
      input_alignment.close();
      return 7;
   }
   input_alignment.close();
   
   //Now that we have the records read in, iterate along the alignment:
   for (size_t i = 0; i < true_one.length(); i++) {
      if (true_one[i] == true_two[i]) { //Homozygous site
         if (test_one[i] == '-' || test_two[i] == '-') { //False indel
            if (test_one[i] != '-') {
               test_one_false_indels++;
               if (position_output_flag) {
                  cout << "False indel at position " << i+1 << endl;
               }
            }
            if (test_two[i] != '-') {
               test_two_false_indels++;
               if (position_output_flag) {
                  cout << "False indel at position " << i+1 << endl;
               }
            }
         } else if (test_one[i] != test_two[i]) {
            if (test_one[i] != true_one[i]) { //False SNP
               test_one_false_snps++;
               if (position_output_flag) {
                  cout << "False SNP at position " << i+1 << endl;
               }
            } else {
               test_two_false_snps++;
               if (position_output_flag) {
                  cout << "False SNP at position " << i+1 << endl;
               }
            }
         }
      } else { //Heterozygous SNP or indel
         if (true_one[i] != '-' && true_two[i] != '-') { //Heterozygous SNP
            //Check the first test haplotype:
            if (test_one[i] == true_one[i]) {
               if (test_one_id == 2) { //Phase switch occurred
                  test_one_switches++;
                  if (position_output_flag) {
                     cout << "Test haplotype 1 switches at position " << i+1 << endl;
                  }
               }
               test_one_id = 1;
            } else if (test_one[i] == true_two[i]) {
               if (test_one_id == 1) { //Phase switch occurred
                  test_one_switches++;
                  if (position_output_flag) {
                     cout << "Test haplotype 1 switches at position " << i+1 << endl;
                  }
               }
               test_one_id = 2;
            } else {
               test_one_bad_calls++;
               if (position_output_flag) {
                  cout << "Test haplotype 1 doesn't match either true haplotype at position " << i+1 << endl;
               }
            }
            //Now check the second test haplotype:
            if (test_two[i] == true_one[i]) {
               if (test_two_id == 2) { //Phase switch occurred
                  test_two_switches++;
                  if (position_output_flag) {
                     cout << "Test haplotype 2 switches at position " << i+1 << endl;
                  }
               }
               test_two_id = 1;
            } else if (test_two[i] == true_two[i]) {
               if (test_two_id == 1) { //Phase switch occurred
                  test_two_switches++;
                  if (position_output_flag) {
                     cout << "Test haplotype 2 switches at position " << i+1 << endl;
                  }
               }
               test_two_id = 2;
            } else {
               test_two_bad_calls++;
               if (position_output_flag) {
                  cout << "Test haplotype 2 doesn't match either true haplotype at position " << i+1 << endl;
               }
            }
         } else { //Indel
            //Not doing anything right now with indels
            if (position_output_flag) {
               cout << "True indel at position " << i+1 << endl;
            }
            if (test_one[i] != true_one[i] && test_one[i] != true_two[i]) {
               test_one_false_snps++;
               if (position_output_flag) {
                  cout << "False SNP due to test haplotype 1 at position " << i+1 << endl;
               }
            } else if (test_two[i] != true_one[i] && test_two[i] != true_two[i]) {
               test_two_false_snps++;
               if (position_output_flag) {
                  cout << "False SNP due to test haplotype 2 at position " << i+1 << endl;
               }
            }
         }
      }
   }
   
   //Output the results:
   cout << "Haplotype switches for test haplotype 1: " << test_one_switches << endl;
   cout << "Haplotype switches for test haplotype 2: " << test_two_switches << endl;
   cout << "False SNPs in haplotype 1: " << test_one_false_snps << endl;
   cout << "False SNPs in haplotype 2: " << test_two_false_snps << endl;
   cout << "False indels in haplotype 1: " << test_one_false_indels << endl;
   cout << "False indels in haplotype 2: " << test_two_false_indels << endl;
   cout << "Bad base calls in haplotype 1: " << test_one_bad_calls << endl;
   cout << "Bad base calls in haplotype 2: " << test_two_bad_calls << endl;
   
   return 0;
}