/*
pagerank weights links based on the amount of links on a page and where they link to
in the wikipedia race, if the starting article does not have a link to the end article, go to a linked article with the highest pagerank
if the desired article is not there, go to the next article with the highest pagerank
so if the previous article had the highest pagerank, return to that page, else continue to the next highest pagerank available
*/

/*
one possible way - power matrix
create a incidence matrix / table representing the graph
create a transition matrix, then recompute it with the proper damping factor
then compute the matrix power to the predetermined convergence factor
well documented - examples galore
*/

/*
Read File
Get amount of verts, create a 2d array with side sizes equal to amount of verts
assign placement based on location in vector
0 if no link going to there, 1 if so
*/
#include "./pagerank.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cstdio>
#include <iterator>


//configure pagerank, run first
pagerank::pagerank() {

    count = 0;
    std::ifstream link_file;
    std::ifstream vertex_file;
    std::string vert;
    std::string edge;

    //read graph from file

    vertex_file.open("data/articles.tsv");
    link_file.open("data/links_without_header.tsv");

    //skip header

    for (int i = 0; i < 12; i++) {
        std::getline(vertex_file, vert);
    }
    //put articles into vector
    while (std::getline(vertex_file, vert)) {
        articles.push_back(vert);
    }
    //put links into a pair then into a vector
    while (std::getline(link_file, edge)) {
        std::string art;
        std::string link;
        
        int ind = edge.find("\t");
        
        art = edge.substr(0, ind);
        
        link = edge.substr(ind + 1);
        
        directed_link = std::pair<std::string, std::string>(art, link);
        links.push_back(directed_link);
    }

    //std::cout << articles.size() << " " << links.size() << std::endl;

    //create adjacency matrix
    adj = new double * [articles.size()];
    for (unsigned i = 0; i < articles.size(); i++) {
        adj[i] = new double[articles.size()];
    }

    //create probability matrix
    transition = new double * [articles.size()];
    for (unsigned i = 0; i < articles.size(); i++) {
        transition[i] = new double[articles.size()];
    }

    //set adjacency matrix equal to 0
    for (unsigned i = 0; i < articles.size(); i++) {
        for (unsigned j = 0; j < articles.size(); j++) {
            adj[i][j] = 0;
        }
    }

    //configure adjacency matrix
    int c = 0;
    int b = 0;
    for (unsigned i = 0; i < articles.size(); i++) {
        
        for (unsigned j = 0; j < links.size(); j++) {
            int index = -1;
            b++;
            //check to see if article has link then set it equal to 1 if so
            if (articles[i] == links[j].first) {
                //std::cout << c++ << std::endl;
                c++;
                auto it = std::find(articles.begin(), articles.end(), links[j].second);
                //std::cout << __LINE__ << std::endl;
                index = it - articles.begin();
                adj[i][index] = 1;
            }
        }
    }    
}

//print out adj
void pagerank::printAdj() {
    if (articles.size() == 0) {
        return;
    }
    for (unsigned i = 0; i < articles.size() - 1; i++) {
        printf( "\n         ");
        for (unsigned j = 0; j < articles.size() - 1; j++) {
            printf( "%1.0f", adj[i][j]);
        }
    }
    printf("\n");
    getchar();
}

//creates transition matrix
void pagerank::createTrans() {
    //if no articles exit
    if (articles.size() == 0) {
        return;
    } 
    
    
    

    // initial transition matrix
    for (unsigned i = 0; i < articles.size() - 1; i++) {
        double sum = 0;
        for (unsigned j = 0; j < articles.size() - 1; j++) sum += adj[i][j];
        if (sum > 0) {
            for (unsigned j = 0; j < articles.size() - 1; j++) {
                //do probabilities
                transition[i][j] = adj[i][j]/sum;
            }
        } else {
            for (unsigned j = 0; j < articles.size() - 1; j++) {
                transition[i][j] = 1.0 / double(articles.size());
            }
        }
    }
    //recomputed transition matrix based on damping factor of .85, generally accepted value
    double alpha = 0.85;
    for (unsigned i = 0; i < articles.size() - 1; i++) {
        for (unsigned j = 0; j < articles.size() - 1; j++) {
            double entry = transition[i][j];
            entry = (alpha * entry) + ((1.0-alpha)/ double(articles.size()));
            transition[i][j] = entry;
        }
    }
}

//prints out transition matrix
void pagerank::printTrans() {
    if (articles.size() == 0) {
        return;
    }
    
     for (unsigned i = 0; i < articles.size() - 1; i++) {
        printf( "\n         ");
        for (unsigned j = 0; j < articles.size() - 1; j++) {
            printf( "%1.4f", transition[i][j]);
        }
    }
    printf("\n\n");
    getchar();
    

}

//finds the pagerank
void pagerank::powerCalc() {
    if (articles.size() == 0) {
        return;
    }
    //std::cout << __LINE__ << std::endl;
    //create the current matrix
    double ** cur = new double*[articles.size()];
    for (unsigned i = 0; i < articles.size(); i++) {
        cur[i] = new double[articles.size()];
    }

    //determine values in current matrix
    for (unsigned i = 0; i < articles.size(); i++) {
        for (unsigned j = 0; j < articles.size(); j++) {
            cur[i][j] = (i == j) ? 1.0 : 0.0;
        } 
    }
    //std::cout << __LINE__ << std::endl;
    int stepC = 0;
    //see if a pagerank vector already exists
    if (pageranks.size() != 0) {
        pageranks.clear();
    }
    //loop through matrices finding the power of the matrix
    do {
        //change max stepC for increased accuracy
        if (stepC >= 4) {
            break;
        }
        
        //create a product matrix
        double ** product = new double*[articles.size()];
        for (unsigned i = 0; i < articles.size(); i++) {
            product[i] = new double[articles.size()];
        }

        //find product of cur and transition matrices
        
        for (unsigned i = 0; i < articles.size(); i++) 
        for (unsigned j = 0; j < articles.size(); j++)
            {
                double summer = 0.0;
                for (unsigned k = 0; k < articles.size(); k++) 
                    summer += cur[i][k] * transition[k][j];
                product[i][j] = summer;
            }

        
        //set current equal to product matrix
        for (unsigned i = 0; i < articles.size(); i++) {
            for (unsigned j = 0; j < articles.size(); j++) {
                cur[i][j] = product[i][j];
            }
        }

        //std::cout << __LINE__ << std::endl;
        //check the square difference of the matrix
        double diff, squarediff = 0.0;
        for (unsigned j = 0; j < articles.size(); j++) {
            for (unsigned i = 1; i < articles.size(); i++) {
                diff = (cur[i][j] - cur[0][j]);
                squarediff += diff * diff;
                
            }
        }
        for (unsigned i = 0; i < articles.size(); i++) {
        delete product[i];
        product[i] = NULL;
        }
        delete[] product;
        //std::cout << squarediff << std::endl;
        //want a very small square diff but takes a very long time due to matrix multiplication of very large matrices
        if (squarediff < 0.00001) break;
        else ++stepC;
    }
    while (1);
    

    std::cout << __LINE__ << std::endl;
    //create pageranks vector
    for (unsigned j = 0; j < articles.size(); j++) {
        pageranks.push_back(cur[0][j]);
    }
    //create pages with ranks vector
    for (unsigned i = 0; i < pageranks.size(); i++) {
        pageWR.push_back(std::pair<std::string, double>(articles[i], pageranks[i]));
    }

    //std::cout << __LINE__ << std::endl;
    //create pageranks file to contain pageranks
    std::ofstream output_file;
    output_file.open("data/pageranks.txt", std::ofstream::out | std::ofstream::trunc);
    std::ostream_iterator<double> output_iterator(output_file, "\n");
    std::copy(pageranks.begin(), pageranks.end(), output_iterator);
}

//find the a path between two articles based on greatest pagerank, not fastest way or best way
int pagerank::pageSearch(std::string startpage, std::string endpage) {

    int startindex, endindex;
    //emergency exit
    if (count > articles.size()) {
        visitedPages.clear();
        return -2;
    }
    //find starting page
    auto it = std::find(articles.begin(), articles.end(), startpage);
    if (it != articles.end()) {

        startindex = it - articles.begin();
    } else {
        visitedPages.clear();
        std::cout << "destination not found" << std::endl;
        return -1;
    }

    //find ending page
    it = std::find(articles.begin(), articles.end(), endpage);
    if (it != articles.end()) {

        endindex = it - articles.begin();
    } else {
        visitedPages.clear();
        std::cout << "destination not found" << std::endl;
        return -1;
    }

    //put starting page into visited pages
    visitedPages.push_back(articles[startindex]);

    //if there is a link between start and end exit
    if (adj[startindex][endindex] == 1) {
        //write visited pages to file
        visitedPages.push_back(articles[endindex]);
        std::ofstream output_file;
        output_file.open("data/searchFile.txt", std::ofstream::out | std::ofstream::trunc);
        std::ostream_iterator<std::string> output_iterator(output_file, "\n");
        std::copy(visitedPages.begin(), visitedPages.end(), output_iterator);
        visitedPages.clear();
        return 0;
    }

    //std::cout << startindex << std::endl;
    //create largeRank which represents the largest rank found
    std::pair<std::string, double> largeRank = std::pair<std::string, double>("", 0);
    
    //loop through articles looking for matching links with the largest pagerank
    for (unsigned i = 0; i < articles.size(); i++) {
        if (adj[startindex][i] == 1) {
            auto it = std::find(visitedPages.begin(), visitedPages.end(), pageWR[i].first);
            if (it == visitedPages.end()) {
                if (pageWR[i].second > largeRank.second) {
                    largeRank = pageWR[i];
                }
            }
        }
    }

    

    count++;
    pageSearch(largeRank.first, endpage);
    return 1;
}

//read pageranks from files
void pagerank::readPageRank() {
    //check if pageranks exists
    if (pageranks.size() != 0) {
        pageranks.clear();
    }
    //check if page with rank exists
    if (pageWR.size() != 0) {
        pageWR.clear();
    }

    std::ifstream pageR;
    std::string rank;
    //read pageranks and put into vectors
    pageR.open("data/pageranks.txt");
    while (std::getline(pageR, rank)) {

        pageranks.push_back(std::stod(rank));
        
    }
    //put into pagewith ranks
    for (unsigned i = 0; i < pageranks.size(); i++) {
        pageWR.push_back(std::pair<std::string, double>(articles[i], pageranks[i]));
    }
}

pagerank::~pagerank() {
    clear();
}

void pagerank::clear() {
    for (unsigned i = 0; i < articles.size(); i++) {
        delete transition[i];
        delete adj[i];
        adj[i] = NULL;
        transition[i] = NULL;
    }
    delete[] transition;
    delete[] adj;
}