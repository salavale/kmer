/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerCounter.cpp
 * @author Juan Salamanca Valenzuela <jsalavale@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "KmerCounter.h"
#include<cmath>
#include<fstream>

using namespace std;

/**
 * DEFAULT_VALID_NUCLEOTIDES is a c-string that contains the set of characters
 * that will be considered as valid nucleotides. 

 * The constructor of the class KmerCounter uses this c-string as a 
 * default parameter. It is possible to use a different c-string if that
 * constructor is used with a different c-string
 */
const char* const KmerCounter::DEFAULT_VALID_NUCLEOTIDES="ACGT";

KmerCounter::KmerCounter(int k, const std::string& validNucleotides):
    _k(k), _validNucleotides(validNucleotides), _allNucleotides(Kmer::MISSING_NUCLEOTIDE + validNucleotides)
{        
    allocate(this->getNumRows(), this->getNumCols());   //We allocate the space we need for the matrix
    initFrequencies();
}

KmerCounter::KmerCounter(const KmerCounter& orig){
    _k=orig._k;
    _validNucleotides = orig._validNucleotides;
    _allNucleotides=orig._allNucleotides;
    allocate(orig.getNumRows(), orig.getNumCols());
    
    for(int i=0; i<orig.getNumRows(); i++){
        for(int j=0; j<orig.getNumCols(); j++){
            _frequency[i][j] = orig._frequency[i][j];
        }
    }
}

KmerCounter::~KmerCounter(){
    deallocate();
    _k = 0;
    _validNucleotides = DEFAULT_VALID_NUCLEOTIDES;
    _allNucleotides = Kmer::MISSING_NUCLEOTIDE + DEFAULT_VALID_NUCLEOTIDES;
}

int KmerCounter::getNumNucleotides() const{
    return _allNucleotides.length();
}

int KmerCounter::getK() const{
    return _k;
}

int KmerCounter::getNumKmers()const{
    return pow(_allNucleotides.length(),_k);
}

int KmerCounter::getNumberActiveKmers()const{
    int active_kmers = 0;
    for(int i=0; i< this->getNumRows();i++){
        for(int j=0; j< this->getNumCols();j++){
            if((*this)(i,j) > 0)
                active_kmers++;
        }
    }
    return active_kmers;
}

std::string KmerCounter::toString() const{
    string outputString = _allNucleotides + " " + to_string(_k) + "\n";
    
    for(int row=0; row< this->getNumRows(); row++){
        for(int col=0; col< this->getNumCols(); col++){
            outputString += std::to_string((*this)(row,col)) + " ";
        }
        outputString += "\n";
    }
    
    return outputString;
}

void KmerCounter::increaseFrequency(const Kmer& kmer, int frequency){
    std::string orig = kmer.toString();
    Kmer aux(kmer);
    aux.normalize(this->_validNucleotides);
    if(orig != aux.toString()) //if they aren't equal, it's because there was at least
                            //one invalid nucleotide. However, it still allows to search
                            //kmers with any missing nucleotides
        throw std::invalid_argument(std::string("void KmerCounter::increaseFrequency(const Kmer& kmer, int frequency):"
                "the given kmer contains any invalid nucleotide.")); 
    
    int row, col;
    this->getRowColumn(kmer, row, col);
    
    _frequency[row][col] += frequency;
}

KmerCounter& KmerCounter::operator=(const KmerCounter& orig){
    if(&orig!=this){
        deallocate();
        _k=orig._k;
        _validNucleotides = orig._validNucleotides;
        _allNucleotides=orig._allNucleotides;
        allocate(orig.getNumRows(), orig.getNumCols());

        for(int i=0; i<orig.getNumRows(); i++){
            for(int j=0; j<orig.getNumCols(); j++){
                _frequency[i][j] = orig._frequency[i][j];
            }
        }
    }
    return *this;
}

KmerCounter& KmerCounter::operator+=(const KmerCounter& kc){
    if(_k != kc._k)
        throw std::invalid_argument(std::string("KmerCounter KmerCounter::operator+=(KmerCounter kc):"
                "kc has a different K (number of nucleotides in kmers)"));
    if(_allNucleotides != kc._allNucleotides)
        throw std::invalid_argument(std::string("KmerCounter KmerCounter::operator+=(KmerCounter kc):"
                "kc has a different set of nucleotides"));
    Kmer kmer; int freq;
    for(int i=0;i<kc.getNumRows();i++){
        for(int j=0; j<kc.getNumCols();j++){
            kmer = kc.getKmer(i,j);
            freq = kc(i,j);
            increaseFrequency(kmer, freq);
        }
    }
    return *this;
}

void KmerCounter::calculateFrequencies(char fileName[]){
    this->initFrequencies();
    ifstream fi;
    fi.open(fileName);
    if(!fi)
        throw std::ios_base::failure(std::string("void KmerCounter::calculateFrequencies(char fileName[]):"
                "could not open fileName."));
    std::string nucleotides_list;
    fi >> nucleotides_list; //Read nucleotides list
    
    this->initFrequencies();
        
    int kmer_used = nucleotides_list.length() -_k + 1;
    for(int i = 0; (i < kmer_used); i++){
        Kmer aux(nucleotides_list.substr(i, _k));
        // Normalize each Kmer in the array
        aux.normalize(this->_validNucleotides);
        
        this->increaseFrequency(aux, 1);
    }
    fi.close();
}

Profile KmerCounter::toProfile(){
    Profile profile(this->getNumberActiveKmers());
    KmerFreq kmerFreq;
    int index = 0;
    
    for(int row=0; row< this->getNumRows(); row++){
        for(int col=0; col< this->getNumCols(); col++){
            if ((*this)(row,col) > 0) {
                kmerFreq.setFrequency((*this)(row,col));
                kmerFreq.setKmer(this->getKmer(row,col));
                profile.at(index)=kmerFreq;
                index++;
            }
        }
    }
    
    
    return profile;
}

int KmerCounter::getNumRows() const{
    return pow(_allNucleotides.length(), ((_k+1)/2));
}

int KmerCounter::getNumCols() const{
    return pow(_allNucleotides.length(), (_k-((_k+1)/2)));
}

int KmerCounter::getIndex(const std::string& kmer) const{
    int index = 0;
    int base = 1;

    for (size_t i = 0; i < kmer.size(); i++) {
        size_t pos = _allNucleotides.find(kmer[kmer.size()-i-1]);
        if (pos == string::npos)
            return -1;
        index += pos * base;
        base *= _allNucleotides.size();
    }
    return index;
}

string KmerCounter::getInvertedIndex(int index, int nCharacters) const {
    string result(nCharacters, Kmer::MISSING_NUCLEOTIDE);

    for (int i = result.size(); i > 0; i--) {
        result[i - 1] = _allNucleotides[index % _allNucleotides.size()];
        index = index / _allNucleotides.size();
    }
    return result;
}

void KmerCounter::getRowColumn(const Kmer& kmer, int& row, int& column){
    std::string kmerLeft, kmerRight;
    int lengthLeft = ((_k+1)/2); //length of the left subkmer
    
    kmerLeft = kmer.toString().substr(0, lengthLeft);  //left subkmer
    kmerRight = kmer.toString().substr(lengthLeft, _k-lengthLeft); //right subkmer
    
    row = this->getIndex(kmerLeft); //row Index
    column = this->getIndex(kmerRight);//col Index
}

Kmer KmerCounter::getKmer(int row, int column) const{
    
    if((row>=this->getNumRows())||(row<0)||(column>=this->getNumCols())||(column<0))
        throw std::invalid_argument(std::string("Kmer KmerCounter::getKmer(int row, int column) const:"
                "row or column is out of bounds."));
    
    string leftKmer, rightKmer;
    leftKmer = getInvertedIndex(row, ((_k+1)/2));
    rightKmer = getInvertedIndex(column, (_k-((_k+1)/2)));
    return Kmer(leftKmer+rightKmer);
}

void KmerCounter::initFrequencies(){
    for(int i=0; i<this->getNumRows(); i++){
        for(int j=0; j<this->getNumCols(); j++){
            _frequency[i][j] = 0;   //We initialize all elements of the matrix to 0;
        }
    }
}

const int& KmerCounter::operator()(int row, int column) const{
    return _frequency[row][column];
}

int& KmerCounter::operator()(int row, int column){
    return _frequency[row][column];
}

//Private Methods

void KmerCounter::allocate(int numRows, int numCols){
    _frequency = new int*[numRows]; //We allocate the number of rows we want
    for(int i=0; i<numRows; i++){
        _frequency[i] = new int[numCols];//each row will have numCols elements
    }
}

void KmerCounter::deallocate(){
    if(_frequency != nullptr){
        for(int i=0; i< this->getNumRows(); i++){
            delete[] _frequency[i];
        }
        delete[] _frequency;
        _frequency = nullptr;
    }
}

void KmerCounter::reallocate(int numRows, int numCols){
    KmerCounter aux = *this;
    
    deallocate();
    allocate(numRows, numCols);
    for(int i=0; i<this->getNumRows();i++)
        for(int j=0;j<this->getNumCols();j++)
            _frequency[i][j] = aux._frequency[i][j];
}
