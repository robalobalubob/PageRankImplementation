#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>



class pagerank {
    public:
    pagerank();
    ~pagerank();
    void clear();
    void printAdj();
    void createTrans();
    void printTrans();
    void powerCalc();
    int pageSearch(std::string startpage, std::string endpage);
    void readPageRank();
    std::vector<double> pageranks;
    std::vector<std::pair<std::string, double>> pageWR;
    std::vector<std::string> noGoodPages;
    std::vector<std::string> articles;
    std::vector<std::string> visitedPages;
    std::vector<std::pair<std::string, std::string>> links;
    std::pair<std::string, std::string> directed_link;
    
    unsigned count;
    
    double ** adj;
    double ** transition;
    
};