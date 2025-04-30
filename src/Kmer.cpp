/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file Kmer.cpp
 * @author Juan Salamanca Valenzuela <jsalavale@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include<iostream>
#include<cctype>
#include<string>
#include "Kmer.h"

using namespace std;

Kmer::Kmer(int k){
    if(k<=0)
        throw std::invalid_argument (std::string("Kmer::Kmer(int k): k <= 0."));
                
    _text = std::string(k, MISSING_NUCLEOTIDE);
}

Kmer::Kmer(const std::string& text){
    if(text.empty())
        throw std::invalid_argument(std::string("Kmer::Kmer(const "
                "std::string& text): text is empty."));
    _text = text;
}

int Kmer::getK() const{
    return _text.size();
}

int Kmer::size() const{
    return _text.size();
}

std::string Kmer::toString() const{
    return _text;
}

const char& Kmer::at(int index) const{
    if((index < 0)||(index >= signed(size())))
        throw std::out_of_range(std::string("const char& Kmer::at(int index) "
                "const : index is not in the range from 0 to k-1 "
                "(both included)."));
    return _text.at(index);
}

char& Kmer::at(int index){
    if((index < 0)||(index >= signed(size())))
        throw std::out_of_range(std::string("char& Kmer::at(int index) "
                ": index is not in the range from 0 to k-1 (both included)."));
    return _text.at(index);
}

void Kmer::toLower(){
    for(int i = 0; i < size(); i++){
           at(i) = tolower(at(i)); 
    }
}

void Kmer::toUpper(){
    for(int i = 0; i < size(); i++){
           at(i) = toupper(at(i)); 
    }
}

void Kmer::normalize(const std::string& validNucleotides){  
    toUpper();
    for(int i = 0; i < size(); i++){
        if(!IsValidNucleotide(at(i), validNucleotides))
            at(i)= MISSING_NUCLEOTIDE;
    }
}

Kmer Kmer::complementary(const std::string& nucleotides, 
        const std::string& complementaryNucleotides) const{
    if(nucleotides.length() != complementaryNucleotides.length())
        throw std::invalid_argument(std::string("Kmer Kmer::complementary("
                "const std::string& nucleotides, const std::string& "
                "complementaryNucleotides) const : nucleotides and "
                "complementaryNucleotides have differents sizes."));
    Kmer complementary_kmer(size());
    for(int i = 0; i < size(); i++){
        if(nucleotides.find(at(i))!= std::string::npos){
            complementary_kmer.at(i) = 
                    complementaryNucleotides.at(nucleotides.find(at(i)));
        }else
            complementary_kmer.at(i) = at(i);
    }
    return complementary_kmer;
}

void Kmer::write(std::ostream& outputStream) const{
    outputStream.write(_text.c_str(), _text.size()+1); //size()+1 will include
                                    //all characters from _text plus '\0' symbol
}

void Kmer::read(std::istream& inputStream){
    _text = "";
    char kmer_char;
    bool cont = true;
    
    inputStream.get(kmer_char);
    if (kmer_char == '\0') {//if kmer_char is '\0', we don´t add it to _text and
                            //finish the search of characters
        cont = false;
    }
    while (inputStream && cont) {
        _text += kmer_char; //if there was no problem getting a character from the
                            //istream and it wasn't '\0', we add it to _text
        
        inputStream.get(kmer_char);
        if (kmer_char == '\0') {//if kmer_char is '\0', we don´t add it to _text and
                                //finish the search of characters
            cont = false;
        }
    }
}

Kmer& Kmer::operator=(const Kmer &orig){
    if(&orig!=this){
        _text = orig._text;
    }
    return *this;
}

bool IsValidNucleotide(char nucleotide, const std::string& validNucleotides){

    return (validNucleotides.find(nucleotide) != std::string::npos);
}

void toLower(Kmer& kmer){
    for(int i = 0; i < kmer.size(); i++){
           kmer.at(i) = tolower(kmer.at(i)); 
    }
}

void ToUpper(Kmer& kmer){
    for(int i = 0; i < kmer.size(); i++){
           kmer.at(i) = toupper(kmer.at(i)); 
    }
}

std::ostream& operator<<(std::ostream& os, const Kmer& kmer){
    os << kmer.toString();
    return os;
}

std::istream& operator>>(std::istream& is, Kmer& kmer){
    std::string text = "";
    is >> text;
    kmer = Kmer(text);
    return is;
}