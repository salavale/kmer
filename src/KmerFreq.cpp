/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Juan Salamanca Valenzuela <jsalavale@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include "KmerFreq.h"

using namespace std;

KmerFreq::KmerFreq():
    _kmer(Kmer()), _frequency(0)
    {   }

const Kmer& KmerFreq::getKmer() const{
    return _kmer;
}

int KmerFreq::getFrequency() const{
    return _frequency;
}

void KmerFreq::setKmer(const Kmer& kmer){
    _kmer = kmer;
}

void KmerFreq::setFrequency(int frequency){
    if(frequency < 0)
        throw std::out_of_range(std::string("KmerFreq::setFrequency(int "
                "frequency): frequency must be > 0"));
    
    _frequency = frequency;
}

std::string KmerFreq::toString() const{
    return(_kmer.toString()+" "+std::to_string(_frequency));
}

void KmerFreq::write(std::ostream& outputStream) const{
    _kmer.write(outputStream);
    outputStream.write(reinterpret_cast<const char*>(&_frequency), sizeof(int));
}

void KmerFreq::read(std::istream& inputStream){
    _kmer.read(inputStream);
    inputStream.read(reinterpret_cast<char*>(&_frequency), sizeof(int));
}

KmerFreq& KmerFreq::operator=(const KmerFreq &orig){
    if(&orig!=this){
        setFrequency(orig._frequency);
        setKmer(orig._kmer);
    }
    return *this;
}

std::ostream& operator<<(std::ostream& os, const KmerFreq& kmerFreq){
    os << kmerFreq.toString();
    return os;
}

std::istream& operator>>(std::istream& is, KmerFreq& kmerFreq){
    Kmer kmer; int freq;
    is >> kmer >> freq;
    kmerFreq.setKmer(kmer); 
    kmerFreq.setFrequency(freq);
    return is;
}

bool operator>(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2){
    bool gt=false;
    if(kmerFreq1.getFrequency()> kmerFreq2.getFrequency()){
        gt=true;
    }
    if(kmerFreq1.getFrequency()==kmerFreq2.getFrequency()){
        if(kmerFreq1.getKmer().toString() < kmerFreq2.getKmer().toString()){ // Preguntar Preguntar Preguntar Preguntar Preguntar Preguntar
            gt=true;
        }
    }
    return gt;
}

bool operator<(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2){
    bool ls = false;
    if(kmerFreq2 > kmerFreq1){
        ls = true;
    }
    return ls;
}

bool operator==(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2){
    bool eq = false;
    bool eqfreq = (kmerFreq1.getFrequency()== kmerFreq2.getFrequency());
    bool eqkmer = (kmerFreq1.getKmer().toString() == kmerFreq2.getKmer().toString());
    if(eqfreq && eqkmer){
        eq=true;
    }
    return eq;
}

bool operator!=(KmerFreq kmerFreq1, KmerFreq kmerFreq2){
    bool neq=false;
    if(!(kmerFreq1==kmerFreq2)){
        neq = true;
    }
    return neq;
}

bool operator<=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2){
    bool leq = false;
    if(kmerFreq1<kmerFreq2){
        leq=true;
    }
    if(kmerFreq1==kmerFreq2){
        leq=true;
    }
    return leq;
}

bool operator>=(const KmerFreq& kmerFreq1, const KmerFreq& kmerFreq2){
    bool geq = false;
    if(kmerFreq1>kmerFreq2){
        geq=true;
    }
    if(kmerFreq1==kmerFreq2){
        geq=true;
    }
    return geq;
}