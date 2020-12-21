// CPCS 324 Algorithms & Data Structures 2
// Reference starter - Basic Flow Network Package
// 2019, Dr. Muhammad Al-Hashimi

// -----------------------------------------------------------------------
// a maximum flow implementation based on simple graph object with 
// minimal linked-list edge implementation and fields
//


var _v = [], _e = []; 
    _v2 = [], _e2 = [];  // note naming conventions in upload guide


// -----------------------------------------------------------------------
function _main()   
{
    // create flow network object
	var GFNet1 = new FNetwork();
    GFNet1.FNet.label = "Figure 10.4 (Levitin, 3rd edition)";

	// read vertices and edges of flow network 
    GFNet1.readNetwork(_v, _e);

    // print graph info
    GFNet1.FNet.printGraph();

    // print before max flow
    GFNet1.printNetwork();

    // apply Edmonds-Karp mehtod on flow network
    GFNet1.edmondsKarp();

    document.write("<br><br>");

    // print after calling Edmonds-Karp mehtod
    GFNet1.printNetwork();
    document.write("<br><br>");

   // print number of augmented paths returned by Edmonds-Karp
    document.write("Num paths: ", GFNet1.numOfPaths, "<br><br>");
    
    //print theoretical max num of augmenting paths for family of instances
    document.write("MAX: ", GFNet1.max, "<br><br>");
    document.write("Average path length (edges): ", GFNet1.avg, "<br><br>");
    
    // ******************************************************************* //
    
    // create flow network object
	var GFNet2 = new FNetwork();
    GFNet2.FNet.label = "Exercise 10.2: 2b (Levitin, 3rd edition)";

	// read vertices and edges of flow network 
    GFNet2.readNetwork(_v2, _e2);

    // print graph info
    GFNet2.FNet.printGraph();

    // print before max flow
    GFNet2.printNetwork();

    // apply Edmonds-Karp mehtod on flow network
    GFNet2.edmondsKarp();

    document.write("<br><br>");

    // print after calling Edmonds-Karp mehtod
    GFNet2.printNetwork();
    document.write("<br><br>");

   // print number of augmented paths returned by Edmonds-Karp
    document.write("Num paths: ", GFNet2.numOfPaths, "<br><br>");
 
    //print theoretical max num of augmenting paths for family of instances
    document.write("MAX: ", GFNet2.max, "<br><br>");
    document.write("Average path length (edges): ", GFNet2.avg, "<br><br>");

        
}


// -----------------------------------------------------------------------
// network object initial requirements: support Edmonds-Karp maximum flow algorithm

function FNetwork()   
{

	// --------------------
    // student property fields next
    
	 /**
	  create flow network graph
    */
   this.FNet = new Graph(); 

   /**
     marking unlabeled vertex with 2 labels 
   */
   this.vertsLabel = [];  

   /**
     flow network is directed
   */
   this.FNet.digraph = true; 

   /**
     flow network is weighted
   */
   this.FNet.weighted = true;

   /**
     Number of paths
   */
   this.numOfPaths = 0;

   /**
     theoretical max number of augmenting paths
   */
   this.max = 0;

   /**
     Average path length
   */
  this.avg = 0;
	
	
	// --------------------
	// student member methods next
	
	// note following are required method names, you are not required to use all of them
	// but must use the name if you choose to have the method

	// accessor methods: getters
	this.edgeFlow=edgeFlowImpl;                // return (get) flow for argument edge i,j 
	this.edgeCap=edgeCapImpl;                  // return capacity for argument edge i,j  
	this.srcVertex=srcVertexImpl;              // return source vertex (or its id, you decide)
	this.sinkVertex=sinkVertexImpl;            // return sink vertex (or its id, you decide)
	this.getFlowVal=getFlowValImpl;            // return current flow *value* in network  // lara
	this.getFlow=getFlowImpl;                  // return current flow as array of {u,v,flow} objects  // lara
	this.inFlow=inFlowImpl;                    // return incoming flow for argument vertex  // lara
	this.outFlow=outFlowImpl;                  // return outgoing flow for argument vertex  // lara
	
	// accessor methods: setters
	this.setEdgeFlow=setEdgeFlowImpl;          // set flow on argument edge (i,j)  // lara
	this.setFlow=setFlowImpl;                  // set flow to argument (including 0) for all edges  // lara
	this.initFlow=initFlowImpl;                // reset flow to 0 for all edges   
	this.setLabel=setLabelImpl;                // set network label (hide Graph code) 
	
	
	// other possibly useful method names
	this.isSrc=isSrcImpl;                      // true if argument is source vertex of network      
	this.isSink=isSinkImpl;                    // true if argument is sink vertex of network   
	this.isEdge=isEdgeImpl;                    // true if argument vertices form an edge ALERT belong to Graph() but leave as test to students  // hadel lara
	this.isBackwardEdge=isBackwardEdgeImpl;    // true if argument vertices form a backward edge // lara
	this.readNetwork=readNetworkImpl;          // input reader method
	this.printNetwork=printNetworkImpl;        // output network including current flow (reference output of final project
	this.edmondsKarp=edmondsKarpImpl;          // implement the Edmonds-Karp algorithm for maximum flow 
	

}

// -----------------------------------------------------------------------
// similar to starter 8
function Vertex(v)
{
    this.label = v.label; // vertex can be labelled
    this.visit = false; //Mark vertex visited or "seen"
    this.adjacent = new List(); //Head pointer of a vertex adjacency list.
	
	this.adjacentByID = adjacentByIdImpl; // Get id of adjacent vertices in an array. 
    this.insertAdjacent = insertAdjacentImpl; // Insert an edge target vertex in adjacency list. An optional edge weight may also be specified.
    this.vertexInfo = vertexInfoImpl; // Get printable vertex info strings
    this.incidentEdge = incidentEdgesImpl; // Get information of incident edges in an array of custom objects  
	
}

// -----------------------------------------------------------------------
// similar to starter 8
function Edge(vert_i,weight)
{
	this.target_v = vert_i;  // Id of edge target vertex 

    //note if-specified idiom add a value, if wegited undefined is nullable (number or null)
	this.weight = !(weight === undefined) ? weight : null;
    this.flow = 0 ;
}


// -----------------------------------------------------------------------
// similar to starter 8 (REMOVE transitive closure/greedy package stuff)
function Graph()
{

	this.verts = []; //Vertex list
	this.ne = 0; // number of edges
    this.nv = 0; // number of vertices
    this.connectedComp = 0; //Number of connected components
    this.weighted = false; //True if graph is weighted, false otherwise 
    this.digraph = false; //True if graph is directed (digraph), false otherwise 
    this.label = ""; //An identification text string tp label graph 
    this.adjMatrix = []; //Graph adjacency matrix, to be created on demand via Graph#makeAdjMatrix
    this.bfs_order = []; //Array of vertex id in BFS order
    this.dfs_push = []; //Array of vertex id in DFS order
    this.dfs_pop = [];  // DFS pop order output array
    
    
	this.readGraph = better_input;   // Graph reader, defaults to internal format
    this.add_edge = addEdgeImpl3; // Insert an edge
    this.bfs = bfsImpl; // Visit connected vertices breadth-first starting at some vertex
    this.componentInfo = componentInfoImpl; //Get printable connectivity info strings
    this.dfs = dfsImpl; // Visit connected vertices depth-first starting at some vertex
    this.isConnected = isConnectedImpl; // Test if graph is connected returning true, otherwise false
    this.listVerts = listVertsImpl; // List graph vertices using info strings returned by Vertex methods
    this.makeAdjMatrix = makeAdjMatrixImpl3; // Create adjacency (or weight, if graph weighted) matrix
    this.makeGraph = makeGraphImpl; //Create a random graph, details depend on implementation
    this.printGraph = printGraphImpl; // Print graph information
    this.topoSearch = topoSearchImpl; //search-like topological traversal of graph vertices
	

}


// -----------------------------------------------------------------------
// functions used by methods of FNetwork and ancillary objects

// -----------------------------------------------------------------------
// begin student code section
// -----------------------------------------------------------------------

// flow network package (REMOVE transitive closure/greedy package stuff)


/**
   implement the Edmonds-Karp algorithm for maximum flow

   @author Hadeel ALoufi
   @implements FNetwork#edmondsKarp
  
*/

function edmondsKarpImpl()
{
    //initialize flow (support method)
    this.initFlow(); 

    //initialize the queue with the source
    var queue = new Queue();
    //srcVertex --> (support method)
    var source = this.srcVertex();

    //label source vertex
    this.setLabel(source,Infinity,"-");
    //add to empty queue
    queue.enqueue(source);

    //not empty queue
    while(!queue.isEmpty())
    {
        //The first vertex id (Front) in the queue
        var i = queue.dequeue();
        
        for(var j=0; j<this.FNet.nv; j++)
        {

            //for every edge from i to j "forward edges"
            if(this.isEdge(i,j))
            {          
                
                //if j is unlabeled then
                if(this.vertsLabel[j] === undefined)
                {
                    //intial residual value
                    var r;

                    //compute residual value --> capacity "u" - flow "x" 
                    //edgeCap, edgeFlow (support method)
                    r = this.edgeCap(i,j) - this.edgeFlow(i,j);
                    

                    //If the residual is positive (we can augment this edge) 
                    if (r>0)
                    {                                          
                        this.setLabel(j, Math.min(this.vertsLabel[i].First_label, r), i+1);
                        queue.enqueue(j);
                    } 
                }
            }        
        }

        //for every edge from j to i "backward edge"
        for(var j=0; j<this.FNet.nv; j++)
        {  
            //isBackwardEdge --> (support method)
            if(this.isBackwardEdge(i,j))
            {
                //if j is unlabeled then
                if(this.vertsLabel[j] === undefined)
                {
                    //If there is a positive flow "x" from j to i 
                    if(this.edgeFlow(j,i)>0)
                    {    
                        var Firstlabel = this.vertsLabel[i].First_label;       
                        this.setLabel(j, Math.min(Firstlabel, this.edgeFlow(j,i)), -1*(i+1));
                        queue.enqueue(j);
                         
                    }
                }
            }        
        }

        var sink = this.sinkVertex();

        //if the sink vertex has been labeled then
        if(!(this.vertsLabel[sink] === undefined))
        {      

            //augment along the augmenting path found

            //increment the number of paths
            this.numOfPaths = this.numOfPaths +1;

            //Start at the sink and move backwards using "second labels "(second labels = sink+1)
            var j = sink+1;

            //The source hasn't been reached
            while(j != 1)
            {                              
                var parent = this.vertsLabel[j-1].Second_label;  
                
                //The second label of vertex j is positive
                if (parent > 0)
                { 
                    flowValue = this.edgeFlow(parent-1 ,j-1) + this.vertsLabel[sink].First_label;
                    this.setEdgeFlow(parent-1 ,j-1 , flowValue);
                }

                //The second label of vertex j is negative
                else
                { 
                    parent = parent*-1;
                    flowValue = this.edgeFlow(j-1 ,parent-1) - this.vertsLabel[sink].First_label;
                    this.setEdgeFlow(j-1, parent-1, flowValue);
                }
            
                j = parent;
            }

            //erase all vertex labels except the ones of the source
            this.vertsLabel = [];
            this.setLabel(source, Infinity, "-");

            //reinitialize the queue with the source
            queue = new Queue();
            queue.enqueue(source);
            
        }
    }
    this.max = (this.FNet.nv*this.FNet.ne)/2;
    this.avg = (this.FNet.ne/this.numOfPaths);
}

/**
   This function sets flow on argument edge by i, j

   @author Lara Waheed

   @implements FNetwork#setEdgeFlow  
   @param {integer} i vertex ID
   @param {integer} j vertex ID
   @param {integer} flowx flow ID
*/
function setEdgeFlowImpl(i,j,flowx)
{
    var v_i = this.FNet.verts[i];
    var adj = v_i.adjacent.traverse();
    for (var i = 0; i < adj.length; i++){
        if (adj[i].target_v == j){
            adj[i].flow = flowx;
        }
    }
}

/**
   set flow to argument (including 0) for all edges

   @author 

   @implements FNetwork#setFlow 
   @param 
*/
function setFlowImpl()
{}


/**
   reset flow to 0 for all edges

   @author Hadeel ALpifi 
   @implements FNetwork#initFlow
*/
function initFlowImpl()
{
    var New_edges;
    for(var k=0; k<this.FNet.nv ; k++){
        New_edges = this.FNet.verts[k].adjacent.traverse();
        for(var j=0; j<New_edges.length; j++){
            New_edges[j].flow = 0;
        }
    }
}


/**
   set network label 

   @author Hadeel ALoufi 
   @implements FNetwork#setLabel 
   @param {integer} vertex vertex id to set its label 
   @param {integer} firstLabel 
   @param {integer} secondLabel 
*/
function setLabelImpl(vertex, firstLabel, secondLabel)
{
    this.vertsLabel[vertex] = {
        First_label: firstLabel,
        Second_label: secondLabel
    }
}

/**
   Return capacity for argument edge i,j

   @author Hadeel ALoufi

   @implements FNetwork#edgeCap
   @param {ineger} i vertex id
   @param {integer} j vertex id
   @return {integer} capacity for argument edge i,j
*/

function edgeCapImpl(i, j)
{
    var x_Vert;
    var x_adjacent;
    var edgeCap;
    x_Vert = this.FNet.verts[i];
    x_adjacent = x_Vert.adjacent.traverse();
    for (var i = 0; i < x_adjacent.length; i++)
    {
         if(x_adjacent[i].target_v == j)
         {
            edgeCap = x_adjacent[i].weight;
            return edgeCap;
         }
    }
}

/**

  Return (get) flow for argument edge i,j

   @author Hadeel ALoufi

   @implements FNetwork#edgeFlow
   @param {ineger} i vertex id
   @param {integer} j vertex id
   @return {integer} flow for argument edge i,j

*/
function edgeFlowImpl(i, j)
{
var x_Vert;
    var x_adjacent;
    var m;
    var edgeFlow;

    x_Vert = this.FNet.verts[i];
    x_adjacent = x_Vert.adjacent.traverse();
    m = x_adjacent.length;

    for (var i = 0; i < m; i++)
    {
         if(x_adjacent[i].target_v == j)
         {
            edgeFlow = x_adjacent[i].flow;
            return edgeFlow;
         }
    }

}

/**
   true if argument is source vertex of network

   @author 
   @implements FNetwork#isSrc 
   @param 
*/
function isSrcImpl()
{}

/**
   true if argument is sink vertex of network

   @author 
   @implements FNetwork#isSink 
   @param 
*/
function isSinkImpl()
{}

/**
   This method returns true if argument vertices
    form an edge ALERT belong to Graph()

   @author Lara Waheed

   @implements FNetwork#isEdge
   @param {ineger} i vertex ID
   @param {integer} j vertex ID
   @return {boolean} true if two vertices are form an edge
*/

function isEdgeImpl(i, j)
{
    return this.FNet.adjMatrix[i][j] != 0;
}

/**
   This method returns true if argument vertices form a backward edge
   @author Lara Waheed
   
   @implements FNetwork#isBackwardEdge
   @param {ineger} i vertex ID
   @param {integer} j vertex ID
   @return {boolean} true if edge is backward 
*/

function isBackwardEdgeImpl(i, j)
{
    return this.FNet.adjMatrix[j][i] != 0;

}


/**
   Input reader method
   @method readNetwork
   @author Ghida'a Bajba 
   @implements FNetwork#readNetwork
   @param {Array} v set of vertices 
   @param {Array} e set of edges
*/

function readNetworkImpl(v, e)
{
    this.FNet.readGraph(v, e);
    this.FNet.makeAdjMatrix();
}

/**
   Output network including current flow 
   @method printNetwork
   @author Ghida'a Bajba 
   @implements FNetwork#printNetwork
*/

function printNetworkImpl() 
{
    for (var i = 0; i < this.FNet.nv; i++)
    {
        var vertex = this.FNet.verts[i];
        var adjacent2 = vertex.adjacent.traverse();
        document.write(i, " {", vertex.label, "}");
        for (var j =0; j <adjacent2.length; j++)
        {
            document.write(" ", i,"-", adjacent2[j].target_v," ", 
            adjacent2[j].weight," ", adjacent2[j].flow);
            document.write( j == adjacent2.length - 1 ? "" : ",");
        }
        document.write("<br>");
    }
}

/**
   return source vertex (or its id, you decide)
   @method srcVertex
   @author Ghida'a Bajba 
   @implements FNetwork#srcVertex
*/
function srcVertexImpl()
{
    return 0;
}

/**
   return sink vertex (or its id, you decide)
   @method sinkVertex
   @author Ghida'a Bajba 
   @implements FNetwork#sinkVertex
*/
function sinkVertexImpl()
{
    return this.FNet.nv -1;
}

/**
   return current flow *value* in network
   @author 
   @implements FNetwork#getFlowVal 
   @param 
*/
function getFlowValImpl()
{}

/**
   return current flow as array of {u,v,flow} objects
   @author 
   @implements FNetwork#getFlow 
   @param 
*/
function getFlowImpl()
{}

/**
   returns incoming flow for argument vertex
   
   @author Lara Waheed
   @implements FNetwork#inFlow 
   @param 
*/
function inFlowImpl()
{}

/**
   return outgoing flow for argument vertex

   @author 
   @implements FNetwork#outFlow 
   @param 
*/
function outFlowImpl()
{}


// -----------------------------------------------------------------------
// published docs section (ref. assignment page)
// similar to starter 8 (strip line comments, leave outline)
// no JSDOC comments in this section
// -----------------------------------------------------------------------
//Add an edge from a vertex pair (original pre Checkpoint 6). 
function addEdgeImpl(u_i, v_i)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u = this.verts[u_i];
    var v = this.verts[v_i];

    // insert (u,v), i.e., insert v (by id) in adjacency list of u
    u.adjacent.insert(v_i);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {
        v.adjacent.insert(u_i);
    }
}

// -----------------------------------------------------------------------
//Add an edge from a vertex pair and optional weight (from Checkpoint 6).
function addEdgeImpl2(u_i, v_i, weight)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u = this.verts[u_i];
    var v = this.verts[v_i];

    // insert (u,v), i.e., insert v in adjacency list of u
    // (first create edge object using v_i as target, then pass object)
    var v_edge = new Edge(v_i);

    if (!(weight === undefined))
    {

        v_edge.weight = weight;
    }

    u.adjacent.insert(v_edge);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {

        var u_edge = new Edge(u_i);

        if (!(weight === undefined))
        {

            u_edge.weight = weight;

        }

        v.adjacent.insert(u_edge);
    }

}

// -----------------------------------------------------------------------
//Add an edge from a vertex pair and an optional weight. A better implementation which relies on a vertex method to handle adjacency details.
function addEdgeImpl3(u_i, v_i, weight)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u = this.verts[u_i];
    var v = this.verts[v_i];

    // insert (u,v), i.e., insert v in adjacency list of u
    // (first create edge object using v_i as target, then pass object)

    u.insertAdjacent(v_i, weight);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)

    if (!this.digraph)
    {

        v.insertAdjacent(u_i, weight);
    }

}

// -----------------------------------------------------------------------
//Output graph properties and list its vertices (replaces better_output).
function printGraphImpl()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "WEIGHTED, " : "", this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ", this.ne, " EDGES:</p>");
    document.write('<br>');
}

// -----------------------------------------------------------------------
//Output a vertex list using descriptive strings returned by the Vertex object 
function listVertsImpl()
{
    var i, v; // local vars
    for (i = 0; i < this.nv; i++)
    {
        v = this.verts[i];
        document.write("VERTEX: ", i, " {", v.label, "} - VISIT: ", v.visit,
            " - ADJACENCY: ", v.adjacentByID(), "<br>");
    }

}

// -----------------------------------------------------------------------
//Input graph data based on default internal format
function better_input(v, e)
{
    // set vertex and edge count fields
    this.nv = v.length;
    this.ne = e.length;

    // input vertices into internal vertex array
    for (var i = 0; i < this.nv; i++)
    {
        this.verts[i] = new Vertex(v[i]);
    }

    //check if the graph is weighted or not
    this.weighted = e[0].w === undefined ? false : true;

    // input vertex pairs from edge list input array
    // remember to pass vertex ids to add_edge()
    for (var i = 0; i < this.ne; i++)
    {
        this.add_edge(e[i].u, e[i].v, e[i].w);
    }

    // double edge count if graph undirected
    if (!this.digraph)
    {
        this.ne = e.length * 2;
    }
}

// -----------------------------------------------------------------------
//Output graph data
function better_output()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "WEIGHTED, " : "", this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ", this.ne, " EDGES:</p>");

    // report connectivity status if available
    document.write(this.connectedComp == 0 ? '<br> no connectivity info' : "DISCONNECTED: ".concat(this.connectedComp));

    // list vertices
    this.listVerts();
}

// -----------------------------------------------------------------------
// Implement a versatile caller for search-like topological traversals of vertices (initially support DFS and BFS).
// Updates and returns the number of connected components Graph#connectedComp. 
//Stores the id of vertices in visit order in Graph#dfs_push or Graph#bfs_order depending on requested search.
function topoSearchImpl(fun)
{
    connectedComp = 0;

    // mark all vertices unvisited
    for (var i = 0; i < this.nv; i++)
    {
        this.verts[i].visit = false;
    }

    // traverse unvisited connected components
    for (var i = 0; i < this.nv; i++)
    {
        if (!this.verts[i].visit)
        {
            fun == 0 ? (this.connectedComp++, this.dfs(i)) : (this.connectedComp++, this.bfs(i));
        }
    }

    return connectedComp;

}

// -----------------------------------------------------------------------
//Recursively traverse a connected component depth-first starting at passed vertex. Inserts visited vertex id in visit order in #dfs_push.
function dfsImpl(v_i)
{
    // get vertex v by its id
    var v = this.verts[v_i];

    // process v 
    v.visit = true;
    this.dfs_push[this.dfs_push.length] = v_i;

    // recursively traverse unvisited adjacent vertices 
    var w = v.adjacentByID();
    for (var i = 0; i < w.length; i++)
    {
        if (!this.verts[w[i]].visit)
        {
            this.dfs(w[i]);
        }
    }

    len_pop = this.dfs_pop.length;
    this.dfs_pop[len_pop] = v_i; 
}

// -----------------------------------------------------------------------
//Traverse a connected component breadth-first starting at passed vertex. Inserts visited vertex id in visit order in #bfs_order.
function bfsImpl(v_i)
{
    // get vertex v by its id
    var v = this.verts[v_i];

    // process v 
    v.visit = true;

    // initialize queue with v
    var q = new Queue();
    q.enqueue(v_i);


    // while queue not empty
    while (!q.isEmpty())
    {

        // dequeue and process a vertex, u
        var u_i = q.dequeue();
        var u = this.verts[u_i];
        this.bfs_order[this.bfs_order.length] = u_i; //fill the bfs_order array when dequeu the vertex

        // queue all unvisited vertices adjacent to u
        var w = u.adjacentByID();
        for (var i = 0; i < w.length; i++)
        {
            if (!this.verts[w[i]].visit)
            {
                this.verts[w[i]].visit = true;
                q.enqueue(w[i]);
            }
        }
    }

}

// -----------------------------------------------------------------------
//Generate adjacency matrix representation of graph - obsolete naïve implementation
function makeAdjMatrixImpl()
{
    // initially the adjacency matrix by zeroz
    for (var i = 0; i < this.nv; i++)
    {

        this.adjMatrix[i] = [];

        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // for each vertex, set 1 for each adjacent  

        var w = this.verts[i].adjacentByID();

        for (var j = 0; j < w.length; j++)
        {
            this.adjMatrix[i][w[j]] = 1;
        }
    }
}

// -----------------------------------------------------------------------
//Generate adjacency matrix representation of graph -  obsolete first weighted implementation
function makeAdjMatrixImpl2()
{
    // initially create row elements and zero the adjacency matrix
    for (var i = 0; i < this.nv; i++)
    {

        this.adjMatrix[i] = [];

        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // for each vertex, set 1 for each adjacent if unweighted, its weight if graph is weighted
        var adj = this.verts[i].adjacent.traverse();
        for (var j = 0; j < adj.length; j++)
        {
            var edge = adj[j];
            this.adjMatrix[i][edge.target_v] = this.weighted ? edge.weight : 1;
        }
    }

}

// -----------------------------------------------------------------------
//Generate an adjacency matrix from internal adjacency lists,This is a better implementation which relies on vertex methods to obtain incident edges information.
function makeAdjMatrixImpl3()
{
    for (var i = 0; i < this.nv; i++)
    {

        this.adjMatrix[i] = [];

        for (var j = 0; j < this.nv; j++)
        {
            this.adjMatrix[i][j] = 0;
        }

        // for each vertex, set 1 for each adjacent if unweighted, its weight if graph is weighted
        var enodes = this.verts[i].incidentEdge();
        for (var j = 0; j < enodes.length; j++)
        {
            var edge_node = enodes[j];
            this.adjMatrix[i][edge_node.adjVert_i] = this.weighted ? edge_node.edgeWeight : 1;
        }
    }
}

// -----------------------------------------------------------------------
//Return connected status based on internal connectivity field 
function isConnectedImpl()
{
    return this.connectedComp == 0 ? true : false;
}

// -----------------------------------------------------------------------
//Report connectivity information if available
function componentInfoImpl()
{
    var checkConnect;
    if (this.connectedComp == 0)
    {
        checkConnect = "no connectivity info <p>";
    }
    else if (this.connectedComp == 1)
    {
        checkConnect = "CONNECTED <p>";
    }
    else
    {
        checkConnect = "DISCONNECTED <p>" + this.connectedComp;
    }
    return checkConnect;
}

// -----------------------------------------------------------------------
//Get id of adjacent vertices in an array
function adjacentByIdImpl()
{
    var traversal = this.adjacent.traverse();

    var traversal_array = [];

    for (var i = 0; i < traversal.length; i++)
    {
        traversal_array[i] = traversal[i].target_v;
    }

    return traversal_array;
}

// -----------------------------------------------------------------------
//Get information of edges incident to a vertex, return array of custom objects containing edge information
function incidentEdgesImpl()
{
    var edges = this.adjacent.traverse(); //get array of incident edges objects

    var edges_info = [];

    for (var i = 0; i < edges.length; i++)
    { //create an array of json objects which contains edge information.
        edges_info[i] = {
            adjVert_i: edges[i].target_v,
            edgeLabel: "",
            edgeWeight: edges[i].weight
        };
    }

    return edges_info;
}

// -----------------------------------------------------------------------
//Get vertex details in a printable string 
function vertexInfoImpl()
{
    var info = "".concat(" {", this.label, "} - VISIT: ", this.visit, " - ADJACENCY: ", this.adjacentByID(), "<br>");

    return info;
}

// -----------------------------------------------------------------------
//Insert a new edge node in the adjacency list of vertex.
function insertAdjacentImpl(v_i, weight)
{
    this.adjacent.insert(new Edge(v_i, weight));
}

// -----------------------------------------------------------------------
//Create a random graph using vertex id as label. If weighted generate random weights between 1 and 50000.
function makeGraphImpl()
{}

// -----------------------------------------------------------------------
