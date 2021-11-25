#include <chrono>
#include <iostream>
#include <stdio.h>
#include <vector>

#include <algorithm>

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#include <sstream>
#include <fstream>
#include <string>

using namespace std;

struct Edge
{
    int src, dest;
    double weight;

    Edge(int src, int dest, double weight) : src(src), dest(dest), weight(weight) { }
};

//the struct for storing final forest
struct Edge_All
{
    int root, src, dest;
    double weight;

    Edge_All(int root, int src, int dest, double weight) : root(root), src(src), dest(dest), weight(weight) { }
};

//used for sorting and counting, at the last step
bool cmp(const Edge_All &a, const Edge_All &b)
{
    return a.root < b.root;
}

//read graph from .mtx file
void read_graph(char* filename, std::vector<Edge> &edge, int &n_rows, int &n_cols) {
    ifstream input(filename);
    if (input.is_open()) {
        std::string line;
        bool skip_line_flag = false;
        while (std::getline(input, line)) {
                                                                   
            // skip lines starting with '%'
            if (line[0] == '%') {
                skip_line_flag = true;
                continue;
            }
            // end when there are no content
    	    if (line.length()==0) {
                    break;
                }

            // skip the line immediately after '%'
            if (skip_line_flag) {
                skip_line_flag = false;
                // get num_vertex and num_edges

                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;
                int counter = 0;
                while ((pos = line.find(delimiter)) != std::string::npos) {
                    token = line.substr(0, pos);
                    if (counter > 2) {
                        counter = 0;
                    }

                    if (counter == 0) {
                        n_rows = stoi(token);
                    }

                    line.erase(0, pos + delimiter.length());
                    counter++;
                }

                n_cols = stoi(line);

                continue;
            }

            std::string delimiter = " ";
            size_t pos = 0;
            std::string token;
            int counter = 0;
            int row, col;
            double val;
            while ((pos = line.find(delimiter)) != std::string::npos) {
                token = line.substr(0, pos);
                if (counter > 2) {
                    counter = 0;
                }

                if (counter == 0) {
                    col = stoi(token);
                }

                if (counter == 1) {
                    row = stoi(token);
                }

                line.erase(0, pos + delimiter.length());
                counter++;
            }
            // std::cout.precision(25);
            // double val = stod(line);
            row = row-1;
            col = col-1;
	    // cout << "val: " << line << endl;
            val = stod(line);
            edge.push_back(Edge(row, col, val));

        }
        input.close();
    }
}

//find the parent node of current node
int find(vector<pair<int,int>>&trees, int i)
{
    // find root and make root as parent of i
    if (trees[i].second != i)
        trees[i].second = find(trees, trees[i].second);
    return trees[i].second;
}

//merge the selected edges
void Union(vector<pair<int,int>>&trees, int a, int b)
{
    int rootA = find(trees, a);
    int rootB = find(trees, b);

    // merge smaller tree to larger one by comparing rank
    if (trees[rootA].first < trees[rootB].first)
        trees[rootA].second = rootB;
    else if (trees[rootA].first > trees[rootB].first)
        trees[rootB].second = rootA;

        // If ranks are same
    else
    {
        trees[rootB].second = rootA;
        trees[rootA].first++;
    }
}

void Boruvkas_function(vector<Edge> edges, int V, int E, vector<double> &edges_lists)
{

    vector<pair<int,int>>trees;
    //cout<<"boruvka begin"<<endl;
    // Create V single-vertex trees
    for (int i = 0; i < V; i++)
    {
        trees.push_back(make_pair(0,i));
    }
    //Initialising two variables
    //TotalTrees stores total no. of trees
    //MST_total_weight stores total weight of MST
    int TotalTrees = V;
    double MST_total_weight = 0.0;


    // cout<<"Edges of MST are :-"<<endl;
    //Loop till only one tree(MST) left
    bool canMerge = true;

    while (canMerge)
    {
        canMerge = false;
        //A vector is created to store smallest edge
        //of each tree. And initialised to -1
        vector<int> smallest_edge(V,-1);

        // Traverse through all edges and update
        // smallest_edge of every tree
        for (int i=0; i<E; i++)
        {
            // Find trees of vertices(s-d) of current edge
            int setA = find(trees, edges[i].src);
            int setB = find(trees, edges[i].dest);
            // cout<<"setA,setB:"<<setA<<","<<setB<<endl;
            // If two vertices of current edge belong to
            //same tree -->continue
            if (setA == setB)
                continue;

                // Else check if current edge is closer to previous
                // smallest_edge edges of setA and setB
            else
            {
                if (smallest_edge[setA] == -1 ||
                    edges[smallest_edge[setA]].weight > edges[i].weight)
                    smallest_edge[setA] = i;

                if (smallest_edge[setB] == -1 ||
                    edges[smallest_edge[setB]].weight > edges[i].weight)
                    smallest_edge[setB] = i;
            }
        }

        //Add edges to MST
        for (int i=0; i<V; i++)
        {
            //if smallest_edge for current set exists
            if (smallest_edge[i] != -1)
            {
                int setA=find(trees, edges[smallest_edge[i]].src);
                int setB=find(trees, edges[smallest_edge[i]].dest);
                //if they belong to same tree -->continue
                if (setA == setB)
                    continue;
                canMerge = true;
                //calculate the total weight of MST
                MST_total_weight += edges[smallest_edge[i]].weight;
                edges_lists.push_back(double(edges[smallest_edge[i]].src));
                edges_lists.push_back(double(edges[smallest_edge[i]].dest));
                edges_lists.push_back(edges[smallest_edge[i]].weight);


                //If two trees are not same then do the union
                //and decrement the no. of trees
                Union(trees, setA, setB);
                TotalTrees--;
    //            cout<<TotalTrees<<endl;
            }
        }
    }
    //Displaying Total weight of  MST
    // cout<<"Total weight of MST is:"<<MST_total_weight<<endl;
}


//same as Boruvkas_function, but this function is for processing the final minimum spanning forest
//edges_all is for storing final data
void Boruvkas_function_all(vector<Edge> edges, int V, int E, vector<double> &edges_lists, vector<Edge_All> &edges_all)
{

    vector<pair<int, int>> trees;    
    for (int i = 0; i < V; i++) {
        trees.push_back(make_pair(0, i));
    }
    
    int TotalTrees = V;
    double MST_total_weight = 0.0;

    bool canMerge = true;

    while (canMerge) {
        canMerge = false;
        vector<int> smallest_edge(V, -1);

        for (int i = 0; i < E; i++) {
            
            int setA = find(trees, edges[i].src);
            int setB = find(trees, edges[i].dest);
            if (setA == setB)
                continue;

            else {
                if (smallest_edge[setA] == -1 ||
                    edges[smallest_edge[setA]].weight > edges[i].weight)
                    smallest_edge[setA] = i;

                if (smallest_edge[setB] == -1 ||
                    edges[smallest_edge[setB]].weight > edges[i].weight)
                    smallest_edge[setB] = i;
            }
        }

        
        for (int i = 0; i < V; i++) {
            
            if (smallest_edge[i] != -1) {
                int setA = find(trees, edges[smallest_edge[i]].src);
                int setB = find(trees, edges[smallest_edge[i]].dest);
                
                if (setA == setB)
                    continue;
                canMerge = true;

                MST_total_weight += edges[smallest_edge[i]].weight;
                edges_lists.push_back(double(edges[smallest_edge[i]].src));
                edges_lists.push_back(double(edges[smallest_edge[i]].dest));
                edges_lists.push_back(edges[smallest_edge[i]].weight);
                
                Union(trees, setA, setB);
                TotalTrees--;
                
            }
        }
    }


    //store the edges information, inculidng parent node index
    for (int i = 0; i < edges_lists.size(); i=i+3) {
        int root = find(trees, edges_lists[i]);
        edges_all.push_back(Edge_All(root, edges_lists[i], edges_lists[i+1], edges_lists[i+2]));
    }
}


int
main(int argc, char **argv)
{
  if (argc != 2)
    {
      fprintf(stderr, "usage: %s <filename>\n", argv[0]);
      return -1;
    }

  
  /* Parallel */
    MPI_Init(&argc, &argv);
    int n_rows, n_cols;
    int numRows, numRows_P;
    int count=0;
    int V = 0; // Number of vertices
    int E = 0; // Number of edges
    std::vector<Edge> edges;
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    read_graph(argv[1], edges, n_rows, n_cols);

    MPI_Barrier(MPI_COMM_WORLD);

    double startTime = MPI_Wtime();
    E = edges.size();
    V = n_rows;
    
    //1st to (P-1)th node process numRows edges
    numRows = ceil(E/world_size)+1;
    //Pth node process numRows_P edges
    numRows_P = E-numRows*(world_size-1);

    std::vector<Edge> edges_p;

    if (world_rank == world_size-1){
        for (int i = 0; i < numRows_P; ++i) {
            count = numRows*world_rank + i;
            edges_p.push_back(Edge(edges[count].src, edges[count].dest, edges[count].weight));
        }

        
        for (int i = 0; i < count; ++i)
        {
            edges_p.push_back(Edge(0, 0, 0.0));
        }

    }
    else{
        for (int i = 0; i < numRows; ++i) {

            count = numRows*world_rank + i;
            edges_p.push_back(Edge(edges[count].src, edges[count].dest, edges[count].weight));
        }

    }
    //edges_p store the sub graph

    MPI_Barrier(MPI_COMM_WORLD);

    vector<double> edges_lists;

    Boruvkas_function(edges_p, V, edges_p.size(), edges_lists);
    //edges_lists store the sub MST

    MPI_Barrier(MPI_COMM_WORLD);

    //---------------for root node to gather data---------------------
        
    int root = 0;
    int edges_lists_size;
    edges_lists_size = edges_lists.size();

    int *cnts;
    cnts = (int *) malloc(world_size * sizeof(int));

    int *offs;
    offs = (int *) malloc(world_size * sizeof(int));
    offs[0] = 0;

    //MPI_Gather for gather each thread has how many edges
    if(world_rank == root){
      
        MPI_Gather(&edges_lists_size, 1, MPI_INT, 
                   cnts, 1, MPI_INT,
                   root, MPI_COMM_WORLD);

        for (int i = 1; i < world_size; i++) { 
            offs[i] = offs[i - 1] + cnts[i - 1]; 
        }
    }
    else{
        MPI_Gather(&edges_lists_size, 1, MPI_INT, 
                   NULL, 0, MPI_INT,
                   root, MPI_COMM_WORLD);
    }

    //MPI_Gatherv for gather edges from each node
    vector<double> all_edges_lists;

    if(world_rank == root){


        int all_edge_size;
        all_edge_size = cnts[world_size-1] + offs[world_size-1];
        all_edges_lists.resize(all_edge_size);

        MPI_Gatherv(&edges_lists.front(), edges_lists_size, MPI_DOUBLE, 
                    &all_edges_lists.front(), cnts, offs, MPI_DOUBLE,
                    root, MPI_COMM_WORLD);


    }
    else{
        MPI_Gatherv(&edges_lists.front(), edges_lists_size, MPI_DOUBLE, 
                    NULL, NULL, NULL, MPI_DOUBLE,
                    root, MPI_COMM_WORLD);
    }
    //------------------------prepare for the final Boruvkas--------------------------

    if (world_rank == root){

        std::vector<Edge> edges_receive;
        for (int i = 0; i < all_edges_lists.size(); i+=3) {
            edges_receive.push_back(Edge(int(all_edges_lists[i]),int(all_edges_lists[i+1]),all_edges_lists[i+2]));
        }

    
        vector<double> edges_lists_final;
        std::vector<Edge_All> edges_all;
	
        Boruvkas_function_all(edges_receive, V, edges_receive.size(), edges_lists_final, edges_all);
	    sort(edges_all.begin(), edges_all.end(), cmp);

	//--------------------count each tree has how many edges--------------------------
        int temp = edges_all[0].root;
        int count = 1;
	vector<int> root_count;
    	for (int i = 1; i < edges_all.size(); ++i) {
            if (edges_all[i].root != temp){
            	root_count.push_back(count);
            	temp = edges_all[i].root;
            	count = 1;
       	      } else{
            	count ++;
       	      }
   	 }
   	 root_count.push_back(count);
	 double endTime = MPI_Wtime();

	//------------------------------print result--------------------------------------
	
         int edge_count = 0;
         for (int i = 0; i < root_count.size(); ++i) {
             cout<<"-------------------------------------------------------------------"<<endl;
             cout<<"number of edges in spanning_tree is "<< root_count[i] <<endl;
             double total_tree_weight = 0.0;
     
             for (int j = 0; j < root_count[i]; ++j) {
                cout<<"source node / destination node / weight :"<<edges_all[edge_count+j].src<<", "<<edges_all[edge_count+j].dest
                << ", " << edges_all[edge_count+j].weight <<endl;
    
                 total_tree_weight = total_tree_weight + edges_all[edge_count+j].weight;
             }
             edge_count = edge_count + root_count[i];
             cout<<"total_tree_weight:"<<total_tree_weight<<endl;
         }
  	 cout<<"-------------------------------------------------------------------"<<endl;
	 cout<<"process time is: "<< endTime - startTime <<endl;
    }

    MPI_Finalize();


  return 0;
}
