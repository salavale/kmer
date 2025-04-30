/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file LEARN.cpp
 * @author Juan Salamanca Valenzuela <jsalavale@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "Profile.h"
#include "KmerCounter.h"
#include <cstdlib>
#include <fstream>

using namespace std;
/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(std::ostream& outputStream) {
    outputStream << "ERROR in LEARN parameters" << endl;
    outputStream << "Run with the following parameters:" << endl;
    outputStream << "LEARN [-t|-b] [-k kValue] [-n nucleotidesSet] [-p profileId] [-o outputFilename] <file1.dna> [<file2.dna> <file3.dna> .... ]" << endl;
    outputStream << endl;
    outputStream << "Parameters:" << endl;
    outputStream << "-t|-b: text mode or binary mode for the output file (-t by default)" << endl;
    outputStream << "-k kValue: number of nucleotides in a kmer (5 by default)" << endl;
    outputStream << "-n nucleotidesSet: set of possible nucleotides in a kmer (ACGT by default). " 
                 << "Note that the characters should be provided in uppercase" << endl;
    outputStream << "-p profileId: profile identifier (unknown by default)" << endl;
    outputStream << "-o outputFilename: name of the output file (output.prf by default)" << endl;
    outputStream << "<file1.dna> <file2.dna> <file3.dna> ....: names of the input files (at least one is mandatory)" << endl;
    outputStream << endl;
    outputStream << "This program learns a profile model from a set of "<< 
            "input DNA files <file1.dna> <file2.dna> <file3.dna> ...." << endl;
    outputStream << endl;
}

/**
 * This program learns a Profile model from a set of input DNA files (file1.dna,
 * file2.dna, ...). The learned Profile object is then zipped (kmers with any 
 * missing nucleotide or with frequency equals to zero will be removed) 
 * and ordered by frequency and saved in 
 * the file outputFilename (or output.prf if the output file is not provided).
 * 
 * Running sintax:
 * > LEARN [-t|-b] [-k kValue] [-n nucleotidesSet] [-p profileId] [-o outputFilename] <file1.dna> [<file2.dna> <file3.dna> ....]
 * 
 * Running example:
 * > LEARN -k 2 -p bug -o /tmp/unknownACGT.prf ../Genomes/unknownACGT.dna
 * 
 * > cat /tmp/unknownACGT.prf
MP-KMER-T-1.0
bug
7
GG 2
AC 1
AG 1
AT 1
CC 1
GA 1
TA 1
 * 
 * @param argc The number of command line parameters
 * @param argv The vector of command line parameters (cstrings)
 * @return 0 If there is no error; a value > 0 if error
 */
int main(int argc, char *argv[]) {   
// Process the main() arguments
    if(argc < 2){
        showEnglishHelp(cerr);
        return 1;
    }
    const char* output = "output.prf";
    int kValue = 5;
    char outputmode = 't';
    std::string profileId = "unknown";
    std::string nucleotidesSet="ACGT";
    std::string cad;
    
    int argument = 1;
    bool continue_reading = false;
    
    if(string(argv[argument]).at(0) == '-'){
        continue_reading = true;
    }
    while(continue_reading){
        switch(argv[argument][1]){
            case 't':
                outputmode='t';
                argument++;
                break;
            case 'b':
                outputmode='b';
                argument++;
                break;
            case 'p':
                cad = string(argv[argument+1]);
                if((cad.size()>=2)&&(cad.at(0)=='\'')&&(cad.at(cad.size()-1)=='\'')){
                    profileId=cad.substr(1,cad.size()-2);     
                }else{
                    profileId=cad;
                }
                argument+=2;
                break;
            case 'k':
                kValue= atoi(argv[argument+1]);
                argument+=2;
                break;
            case 'n':
                nucleotidesSet=string(argv[argument+1]);
                argument+=2;
                break;
            case 'o':
                output=argv[argument+1];
                argument+=2;
                break;
            default:
                showEnglishHelp(cerr);
                return 1;
        }
        if(argument<argc){   
            if(string(argv[argument]).at(0) != '-'){
                continue_reading=false;
            }
        }else{
            showEnglishHelp(cerr);
            return 1;
        }
    }    
    // Loop to calculate the kmer frecuencies of the input genome files using 
    // a KmerCounter object
    int num_files=argc-argument;
    
    KmerCounter input(kValue, nucleotidesSet);
    KmerCounter total(kValue, nucleotidesSet);
    for(int i=0; i < num_files; i++){
        input.calculateFrequencies(argv[argument+i]);
        total+=input;
    }
    
//    cout << total.toString();
    // Obtain a Profile object from the KmerCounter object
    Profile prf=total.toProfile();
//    cout << prf.toString();
    prf.setProfileId(profileId);
    
    // Zip the Profile object
    prf.zip(true);
    // Sort the Profile object
    prf.sort();
    // Save the Profile object in the output file
    prf.save(output, outputmode);
    
    return 0;
}

