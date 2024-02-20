#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <regex>


// samtools mpileup
int mpileF(std::string name, std::string inpath) {
    std::string command = "samtools mpileup " + inpath + name + ".bam > temp/test.mpileup --no-BAQ --ff UNMAP,DUP -A -Q 0 -q 0";

    int result = system(command.c_str());

    if (result == 0) {
        Rcpp::Rcout  << "mpileup executed successfully." << "\n";
    } else {
        Rcpp::Rcout << "mpileup execution failed." << "\n";
    }

    return 0;
}


// Extract fifth field from mpileup output
std::string extractFifthFieldF(const std::string& line) {

    std::istringstream iss(line);
    std::string field;

    // Skip the first four fields
    for (int i = 0; i < 3; ++i) {
      std::getline(iss, field, '\t');
    }

    // Extract the fifth field
    std::getline(iss, field, '\t');

    return field;
}

// mainBash - count ACGT, compare sums, remove indels and deletions
int mainBashF(std::string name, std::string outpath) {

    // Count number of ACGTs
    std::ifstream inputFile1("temp/test.mpileup");
    std::ofstream cleanFile("temp/clean.mpileup");

    if (!inputFile1.is_open() || !cleanFile.is_open() ) {
      Rcpp::Rcerr << "Failed to open files." << "\n";
      return 1;
    }

    std::string line1;
    while (std::getline(inputFile1, line1)) {
      std::string nucleotides = extractFifthFieldF(line1.substr(line1.find('\t')+1));
      // Count deletions and insertions
      std::regex patterndelt("[^-]");
      std::string delOnly = std::regex_replace(nucleotides, patterndelt, "");

      std::regex patterninst("[^+]");
      std::string instOnly = std::regex_replace(nucleotides, patterninst, "");
      double mount = 0.0;
      mount += (delOnly.size() + instOnly.size());
      //Rcpp::Rcout << mount << std::endl;
      //std::string mout;
      //std::getline(inputFile1, mout);

      if(mount == 0.0){
        cleanFile << line1 << "\n";
      } else{
      }
    }

    inputFile1.close();
    cleanFile.close();

    // Count number of ACGTs
    std::ifstream inputFile("temp/clean.mpileup");
    std::ofstream outputCounts("temp/temp_counts.txt");

    if (!inputFile.is_open() || !outputCounts.is_open() ) {
      Rcpp::Rcerr << "Failed to open files." << "\n";
      return 1;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
      std::string nucleotides = extractFifthFieldF(line.substr(line.find('\t')+1));

      // Count A
      std::regex patterna("[^Aa]");
      std::string aOnly = std::regex_replace(nucleotides, patterna, "");

      // Count C
      std::regex patternc("[^Cc]");
      std::string cOnly = std::regex_replace(nucleotides, patternc, "");

      // Count G
      std::regex patterng("[^Gg]");
      std::string gOnly = std::regex_replace(nucleotides, patterng, "");

      // Count T
      std::regex patternt("[^Tt]");
      std::string tOnly = std::regex_replace(nucleotides, patternt, "");

      // Make sum of counts
      double total = (aOnly.size() + cOnly.size() + gOnly.size() + tOnly.size());

      // Extract columns 1, 2, and 4 from mpileup output
      // This is the contig/chromosome, position, and depth/coverage
      std::istringstream iss(line);
      std::vector<std::string> columns;
      std::string field;
      for (int i = 0; i < 4; ++i) {
        std::getline(iss, field, '\t');
        if (i == 0 || i == 1 || i == 3) {
          columns.push_back(field);
        }
      }

      // Write the extracted columns to the output file
      for (const auto& column : columns) {
        outputCounts << column << '\t';
      }
      outputCounts <<  aOnly.size() << '\t'<< cOnly.size() << '\t' << gOnly.size() << '\t' << tOnly.size() << '\t' << total << '\n';

    }


    inputFile.close();
    outputCounts.close();


    // Compare sum of ACGTs to mpileup & combine
    // Remove INDELS & potential errors in calculation (mpileup != ACGTs sum)

      std::ifstream inputCounts("temp/temp_counts.txt");
      std::ofstream outputCompare("temp/temp_compare.txt");
      std::ofstream outputCounts4("temp/temp_counts4.txt");

       if (!inputCounts.is_open() || !outputCompare.is_open() || !outputCounts4.is_open()) {
                    Rcpp::Rcerr << "Failed to open files." << "\n";
        return 1;
      }


       std::string line2;
       while (std::getline(inputCounts, line2)) {

         std::istringstream iss(line2);
        // std::string out;
         //std::getline(inputCounts, out);


         std::string field;
         int column3;
         int column8;

         // Extract values from columns 3 and 8
         for (int i = 0; i < 8; ++i) {
           if (i == 2) {
             iss >> column3;
           } else if (i == 7) {
             iss >> column8;
           } else {
             iss >> field;
           }
         }

         // Compare values and write the result to the output file
         if (column3 == column8) {
           outputCompare << "TRUE" << '\n';
           outputCounts4 << line2 << "\n";
         } else {
           outputCompare << "FALSE" << '\n';
         }

       }

    inputCounts.close();
    outputCompare.close();
    outputCounts4.close();

    // Extract relevant columns and add header
     std::ifstream inputCounts4("temp/temp_counts4.txt");
     std::ofstream outputFinal(""+ outpath + name + ".txt");

    if (!inputCounts4.is_open() || !outputFinal.is_open()) {
       Rcpp::Rcerr << "Failed to open files." << "\n";
       return 1;
     }

     outputFinal << "chr\tpos\tdepth\tA\tC\tG\tT\n";

      std::string countsLine;
    while (std::getline(inputCounts4, countsLine)) {

      std::istringstream iss2(countsLine);
      std::vector<std::string> columns2;
      std::string field2;
      for (int i = 0; i < 7; ++i) {
        std::getline(iss2, field2, '\t');
        if (i == 0 || i == 1 || i ==2 || i == 3 || i == 4 || i == 5 || i == 6) {
          columns2.push_back(field2);
        }
      }

      // Write the extracted columns to the output file
      for (const auto& column2 : columns2) {
        outputFinal << column2 << '\t';
      }
      outputFinal << '\n';
    }

     inputCounts4.close();
     outputFinal.close();


    // Remove Temporary Files
       //std::remove("temp/temp_counts.txt");
      // std::remove("temp/temp_compare.txt");
       //std::remove("temp/temp_counts4.txt");
       //std::remove("temp/clean.mpileup");

    Rcpp::Rcout << "mpileup is now readable!." << "\n";

    return 0;
 }

//' @title Prepare data - Step 1
//'
//' @description This function transforms a BAM file into a text file.
//'   Specifically, this function uses [samtools mpileup](http://www.htslib.org/doc/samtools-mpileup.html)
//'   to translate your BAM into a tab-separated file. We then filter this file to remove indels and deletions.
//'   When running this function, a temporary folder will be created (named 'temp/'), however this folder will be
//'   removed once the process is complete.
//'
//' @details Warning, due to the processing time needed for samtools mpileup,
//'   this step may take some time. This function also requires samtools to be located locally. Please see
//'   our [Data Preparation](mlgaynor.com/nQuack/DataPreparation.html) article for more information. Warning, this writes a temporary folder
//'   titled 'temp'. If you want to run multiple samples at once, we suggest you set the working directory to separate locations to ensure that
//'   your temp folder/files are not overwritten.
//'
//' @param name File name without the suffix. For example, if your file is called "frog.bam", this input should be "frog".
//' @param inpath Location of input file.
//' @param outpath Location for output file.
//'
//' @returns Writes text file with the following columns: chromosome, position, depth, A, C, G, and T.
//'
// [[Rcpp::export]]
void prepare_data(std::string name, std::string inpath, std::string outpath) {
   std::string command1 = "mkdir temp";
   int result1 =  system(command1.c_str());
   mpileF(name, inpath);
   mainBashF(name, outpath);
   std::string command2 = "rm -rf temp";
   int result2 = system(command2.c_str());
   std::string command3 ="perl -pi -e 's/\t\n/\n/g' -i " + outpath + name + ".txt";
   int result3 = system(command3.c_str());
   if(result1 + result2+ result3 == 1){
      Rcpp::Rcout  << "done" << "\n";
   }

}
