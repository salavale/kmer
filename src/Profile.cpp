/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Juan Salamanca Valenzuela <jsalavale@correo.ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */
#include "Profile.h"
#include <fstream>
#include <cmath>

using namespace std;

const string Profile::MAGIC_STRING_T="MP-KMER-T-1.0";
const string Profile::MAGIC_STRING_B="MP-KMER-B-1.0";

Profile::Profile(){
    _profileId = "unknown";
    _capacity = INITIAL_CAPACITY;
    _size=0;
    allocate(_capacity);
}

Profile::Profile(int size){
    if(size<0)
        throw std::out_of_range(string("Profile::Profile(int size): size must be"
            " greater than 0 "));
    _capacity = size;
    _size = size;
    allocate(_size);
    _profileId = "unknown";
}

Profile::Profile(const Profile& orig){
     _capacity = orig._capacity;
    _size = orig._size;
    allocate(_size);
    _profileId = orig._profileId;
    
    for(int i=0; i<_size;i++){
        _vectorKmerFreq[i]=orig._vectorKmerFreq[i];
    }
}

Profile::~Profile(){
    deallocate();
    _size=0;
    _capacity = INITIAL_CAPACITY;
    
}

Profile& Profile::operator=(const Profile &orig){
    if(&orig!=this){
        deallocate();
        _capacity = orig._capacity;
        allocate(_capacity);
        _size=orig._size;
        _profileId = orig._profileId;
        
        for(int i=0; i<_size; ++i){
            _vectorKmerFreq[i] = orig._vectorKmerFreq[i];
        }
    }
    return *this;
}

const std::string& Profile::getProfileId() const{
    return _profileId;
}

void Profile::setProfileId(const std::string& id){
    _profileId = id;
}

const KmerFreq& Profile::at(int index) const{
    if((index<0)||(index>=_size))
        throw std::out_of_range(std::string("const KmerFreq& Profile::at(int index) const: index must be between"
                " 0 and _size-1"));
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::at(int index){
    if((index<0)||(index>=_size))
        throw std::out_of_range(std::string("KmerFreq& Profile::at(int index): index must be between "
                "0 and _size-1"));
    return _vectorKmerFreq[index];
}

int Profile::getSize() const{
    return _size;
}

int Profile::getCapacity() const{
    return _capacity;
}

const double Profile::getDistance(const Profile& otherProfile) const{
    
    if((_size==0)||(otherProfile.getSize()==0))
        throw std::invalid_argument(std::string("const double Profile::getDistance(const Profile& otherProfile) const: "
                "One of the profiles given to getDistance() was empty"));
    double distance = 0.0;
    
    for(int i=0; i<_size;i++){
        //we search for the kmer at i-rank on the other prf
        int pos_other = otherProfile.findKmer(at(i).getKmer());
        if(pos_other != -1){ //if !=-1 then we have found it
            distance += abs(i - pos_other);
        }else{//if not, we use the size of the other prf
            distance += abs(i - otherProfile.getSize());
        }
    }
  
    return (distance/(_size*otherProfile.getSize()));
}

int Profile::findKmer(const Kmer& kmer, int initialPos, int finalPos)const{
    
    bool not_found = true;
    int pos = initialPos;
    
    if(initialPos<=finalPos){
        while(not_found&&(pos <= finalPos)){
            if(_vectorKmerFreq[pos].getKmer().toString() == kmer.toString())
                not_found = false;
            else
                pos++;
        }
    }
    
    if(not_found)
        pos = -1;
    
    pos = pos-initialPos;
    return pos;
}

int Profile::findKmer(const Kmer& kmer)const{
    return (findKmer(kmer, 0, _size-1));
}

std::string Profile::toString() const{
    std::string str;
    str =  _profileId + '\n';
    str += to_string(_size)+'\n';
    for(int i =0; i<_size;i++){
        str +=  _vectorKmerFreq[i].toString()+'\n'; 
    }
    return str;
}

void Profile::sort(){
    //using insertion sort method
    
    //we start to sort the array from the second element because the first one
    //will be in it place since there is no other element before itself
    for (int left = 1; left < _size; left++) {
        //we saved the data from the element we want to sort
        KmerFreq hold = _vectorKmerFreq[left];
        
        //We search the place of "sorting" and "hold" by moving elements of the 
        //array to the right if sorting is bigger than the elements freq
        int i = left;   //we start the search from the place where was
                        //the element we saved and we go looking the previous elements

        while (i > 0 && (hold > _vectorKmerFreq[i-1])){
            _vectorKmerFreq[i]=_vectorKmerFreq[i-1];
            i--;
        }
        
        //when we finish the search, we insert the KmerFreq we saved before in the
        //place it belong in the sorted array at that cycle
        _vectorKmerFreq[i]=hold;    
        
    } // for
    //at the end of the for cycle, we will have had all the elements at the array
    //gone through the sorting process
}

void Profile::save(const char fileName[], char mode) const{ //Pendiente modificar
    if((mode != 't') && (mode != 'b')){
        throw std::invalid_argument(std::string("void Profile::save(char* fileName[]) const:"
                "The given mode \"" + std::string(1,mode) +
                "\" must be 't' for text mode or 'b' for binary mode."));
    }
    ofstream write(fileName, ios::out|ios::binary);
    if(!write)
        throw std::ios_base::failure(std::string("void Profile::save(char* fileName[]) const:"
                " it wans't possible to open current file."));
    if(mode == 't'){
        write << MAGIC_STRING_T << "\n";
        write << (*this);
        if(!write){
            write.close();
            throw std::ios_base::failure(std::string("void Profile::save(char* fileName[]) const:"
                    " it couldn't write on the file."));
        }
    }else if(mode == 'b'){
        write << MAGIC_STRING_B << "\n";
        write << this->getProfileId() << "\n";
        write << std::to_string(this->getSize()) << "\n";
        for(int i=0 ;i<getSize(); i++){
            _vectorKmerFreq[i].write(write);
        }
    }
    write.close();
}

void Profile::load(const char fileName[]){
    ifstream stream(fileName, ios::in|ios::binary);
    if(!stream){
        throw std::ios_base::failure(std::string("void Profile::load(const char fileName[]):"
            "couldn't open the file"));
    }
    
    std::string str;
    //reads the magic string
    stream >> str;
    stream.ignore();
    if((str != MAGIC_STRING_T)&&(str != MAGIC_STRING_B)){
        stream.close();
        throw std::invalid_argument(std::string("void Profile::load(const char fileName[]): invalid magic string was found in the given file"));
    }
    if(str == MAGIC_STRING_T){
        stream >> *this;
    }else if(str == MAGIC_STRING_B){
        //reads the profileId
        std::string prf_id;
        getline(stream, prf_id);
        //Read the number of kmers that are on the file
        int nkmers;
        stream >> nkmers;
        if(nkmers<0){
            stream.close();
            throw std::out_of_range(std::string("void Profile::load(const char fileName[]): "
                    "nkmers must be greater than 0"));
        }
        
        *this = Profile(nkmers);
        this->setProfileId(prf_id);
        
        KmerFreq aux;
        for(int i = 0; i < nkmers; i++){
            aux.read(stream);
            _vectorKmerFreq[i] = aux;
        }
    }
    if(!stream)
    {
        stream.close();
        throw std::ios_base::failure(std::string("void Profile::load(const char fileName[]):  an error occurs while reading from the file"));
    }
    stream.close();
}

void Profile::append(const KmerFreq& kmerFreq){
    int pos = findKmer(kmerFreq.getKmer());
    if(pos == -1){
        if(_size==_capacity){
            
        reallocate(_capacity+BLOCK_SIZE);   
        }
        _vectorKmerFreq[_size] = kmerFreq;
        _size++;
    }else{
        _vectorKmerFreq[pos].setFrequency(_vectorKmerFreq[pos].getFrequency()+kmerFreq.getFrequency());
    }
}

void Profile::normalize(const std::string& validNucleotides){
    // Loop to traverse and normalize each one of the kmers in array
    for(int i=0; i<_size; i++){
        Kmer k = _vectorKmerFreq[i].getKmer();
        k.normalize(validNucleotides); // Normalize kmer i
        _vectorKmerFreq[i].setKmer(k);
    }
    // Loop to traverse the kmers in array from position 1 to position size-1
    int i = 1;
    while(i<_size){
        int pos = findKmer(_vectorKmerFreq[i].getKmer(), 0, i-1);
        if(pos != -1){
            _vectorKmerFreq[pos].setFrequency(_vectorKmerFreq[pos].getFrequency()+_vectorKmerFreq[i].getFrequency());
            deletePos(i);
        }else{
            i++;
        }
    }
}

void Profile::deletePos(int pos){
    if((0>pos)||(pos>=_size))
        throw std::out_of_range(std::string("DeletePosArrayKmerFreq(KmerFreq "
                "array[], int &nElements, int pos)pos("+std::to_string(pos)+")"
                " must be between 0 and nElements-1("+std::to_string(_size-1)+")"));
    for(int i=pos; i<_size-1; i++){
        _vectorKmerFreq[i] = _vectorKmerFreq[i+1];
    }    
    _size--;
}

void Profile::zip(bool deleteMissing, int lowerBound){
    int i = 0;
    while (i < _size){
        if (((_vectorKmerFreq[i].getKmer().toString().find(Kmer::MISSING_NUCLEOTIDE) != std::string::npos)&&deleteMissing)
                ||(_vectorKmerFreq[i].getFrequency() <= lowerBound)){
            deletePos(i);
        } else {
            i++;
        }
    }
}

void Profile::join(const Profile& profile){
    for(int i=0; i<profile.getSize(); i++){
        append(profile.at(i));
    }
}

KmerFreq& Profile::operator[](int index) const{ //No Seguro No Seguro No Seguro No Seguro
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::operator[](int index){ // No Seguro No Seguro No Seguro No Seguro No Seguro
    return _vectorKmerFreq[index];
}


Profile& Profile::operator+=(const KmerFreq& kmerFreq){
    append(kmerFreq);
    return *this;    // No seguro No seguro No seguro No seguro No seguro
}

Profile& Profile::operator+=(const Profile& profile){
    for(int i=0; i<profile.getSize(); i++){
        append(profile[i]);
    }
    return *this; 
}

std::ostream& operator<<(std::ostream& os, const Profile& profile){
    os << profile.toString();
    return os;
}

std::istream& operator>>(std::istream& is, Profile& profile){
    
    std::string id;
    getline(is,id);
    int nkmerfreqs;
    is >> nkmerfreqs;
    
    profile = Profile(nkmerfreqs);
    profile.setProfileId(id);

    KmerFreq aux;
    for(int i=0; i<nkmerfreqs; i++){
        //is >> aux;
        is >> aux;
        profile.at(i) = aux;
    }
    return is;
}

//Private Methods

void Profile::allocate(int nelements){
    _vectorKmerFreq = new KmerFreq[nelements];
}

void Profile::deallocate(){
    if(_vectorKmerFreq != nullptr){
        delete[] _vectorKmerFreq;
        _vectorKmerFreq = nullptr;
    }
}

void Profile::reallocate(int nelements){
    Profile aux = *this;
    
    deallocate();
    allocate(nelements);
    for(int i=0; i<aux._size;i++)
        _vectorKmerFreq[i] = aux._vectorKmerFreq[i];
}