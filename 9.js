// CPCS 324 Algorithms & Data Structures 2
// Graph data structure starter - Transitive Closure Package
// 2020, Dr. Muhammad Al-Hashimi

// -----------------------------------------------------------------------
// simple graph object with linked-list edge implementation and minimal fields
// extra vertex and edge member fields and methods to be added later as needed
//

var _v = [], _e = [];   


// -----------------------------------------------------------------------

function _main()   
{
    
	// create an undirected graph 

    var g = new Graph();
    

    // set input graph properties (label, directed etc.)

    g.label = 'Exercise 9.2: 1b (Levitin, 3rd edition)';


    // use global input arrays _v and _e to initialize its internal data structures

    g.readGraph(_v, _e);

    g.makeAdjMatrix(); 


    // use printGraph() method to check graph

    g.printGraph();


    // perform depth-first search and output stored result

    g.connectedComp = g.topoSearch(g.dfs);

    document.write("dfs_push: ", g.dfs_push, "<br><br>");


    // report connectivity status if available

    document.write(g.componentInfo());


    // perform breadth-first search and output stored result

    g.connectedComp = g.topoSearch(g.bfs);

    document.write("bfs_order: ", g.bfs_order, "<br>");


    
	// output transitive closure matrix using Warshall-Floyd
    
    document.write('<br>Transitive closure<br>');
    
    g.warshallFloyd();
    
    if(g.adjMatrix[0][0].tc === undefined)
    {
        document.write("Transitive closure is computed for directed graph only<br>");
    }
    else
    {
        for(var i=0; i<g.nv; i++){
            for(var j=0; j<g.nv; j++){
                document.write(g.adjMatrix[i][j].tc,',');
            }
            document.write('<br>');
        }
    }


    // output distance matrix using Warshall-Floyd
	
	document.write('<br>Distance matrix <br>');
    
    g.warshallFloyd();
    	
    if(g.adjMatrix[0][0].dist === undefined)
    {
        document.write("Distance matrix is computed for weighted graph only<br>");
    }
    else
    {
        for(var i=0; i<g.nv; i++){
            for(var j=0; j<g.nv; j++){
                document.write(g.adjMatrix[i][j].dist,',');
            }
            document.write('<br>');
        }
    }

    
    //print minimal spanning tree

    document.write('<br>MST by Prim2 (linear PQ)<br>');
    
    g.prim();
    if(g.weighted){
        printTree(g.VT);
    }
    else
    {
        document.write('minimal spanning tree works for weighted graph only<br>');
    }

    //Print shortest path by Dijkstra from vertex 0. Shortest path problems works for weighted graph only

    document.write('<br>Shortest paths by Dijkstra from vertex 0<br>');
    g.Dijkstra(0);
    if(g.weighted){
        printTree(g.VT);
    }
    else
    {
        document.write('Shortest path problems works for weighted graph only<br>');
    }

    
    //Print distance matrix from Dijkstra. Distance matrix is computed for weighted graph only
    document.write('<br>Distance matrix from Dijkstra<br>');
    
    var disMDijkstra = [];
    
    if(g.weighted){
        for(i=0; i<g.nv; i++){
            //The shortes path from vertex i to all other vertices
            g.Dijkstra(i);             
            disMDijkstra[i] = [];

            //fill the distances of row i
            for(j=0; j<g.nv; j++){
                disMDijkstra[i][g.VT[j].tree] = g.VT[j].distance;
            }

            //print row i
            document.write(disMDijkstra[i],'<br>');
        }
    }
    else
    {
        document.write('Distance matrix is computed for weighted graph only<br>');
    }
}
// -----------------------------------------------------------------------

function Vertex(v)
{
	// published docs section (ref. specs, assignment page)
	// for this section, strip line comments
	// no JSDOC comments in this section
	
	// base property fields from P1M1
	
	// base property fields
    this.label = v.label; // vertex can be labelled
    this.visit = false; // vertex can be marked visited or "seen"
    this.adjacent = new List(); // init an adjacency list

    // base member methods
    this.adjacentByID = adjacentByIDImpl; // Get id of adjacent vertices in an array.
    this.incidentEdges = incidentEdgesImpl; // return target id of incident edges in array
    this.vertexInfo = vertexInfoImpl; // Get vertex details in a printable string
    this.insertAdjacent = insertAdjacentImpl; // Insert a new edge node in the adjacency list of vertex.	

		
	// -------------------------------------------
	// more student fields next
	
	
	// --------------------
	// more student methods next
	
}

// -----------------------------------------------------------------------

function Edge(vert_i,weight)
{
	// published docs section (ref. specs, assignment page)
	// for this section, strip (delete) line comments
	// no JSDOC comments in this section
   // --------------------
	// base property fields
    this.target_v = vert_i; //Id of edge target vertex
    this.weight = weight; //Edge weight/cost
    
    
	
	// base member methods

	// -------------------------------------------
	// more student fields next
	
	
	// --------------------
	// more student methods next

}


// -----------------------------------------------------------------------

function Graph()
{
	// published docs section (ref. specs, assignment page)
	// for this section, strip (delete) line comments
	// no JSDOC comments in this section
   // --------------------
	// base property fields

	this.verts = [];
	this.nv = 0;  // ... etc. from P1M1 (remove)
	this.ne = 0;
    this.digraph = false;
    this.weighted = false;
    this.dfs_push = [];
    this.bfs_out = [];
    this.label = "";
    this.nc = 0;
    this.adjMatrix = [];
    this.VT = []; 
	

	// base member methods
	
	this.readGraph = better_input;  // ... (complete next)
	this.add_edge = addEdgeImpl2;
	this.printGraph = printGraphImpl;
	this.makeGraph = makeGraphImpl;
	this.list_verts =listVertsImpl;
	this.dfs = dfsImpl;
	this.bfs = bfsImpl;
	this.makeAdjMatrix = makeAdjMatrixImpl2;
	this.isConnected = isConnectedImpl;
	this.componentInfo = componentInfoImpl;
	this.topoSearch = topoSearchImpl;
	
	
	// -------------------------------------------
	// more student fields next

    // --------------------
     // transitive closure package 
    /**
    	Distance matrix of shortest paths, set after applying warshallFloyd.
    	Stores output of warshallFloyd call. 
    	@default [ ]
    */
   this.floydD = [];
   /**
       Transitive closure matrix , set after applying warshallFloyd.
       Stores output of last warshallFloyd call.
       @default [ ]
   */
   this.warshallTC = [];
   /**
       Transitive closure matrix , set after applying DfsTC.
       Stores output of last DfsTC call.
       @default [ ]
   */
   this.dfsTC = [];
   /**
       MST list of vertices , set after applying PrimMST.
       Stores output of last PrimMST call.
       @default [ ]
   */
   this.Prim_Edge = [];

   // student methods
   /**
       Check if there is a path between two vertices in a digraph
       @method
   */
   this.hasPath = hasPathImpl;
   /**
       Get the length of the shortest path between two vertices in a weighted graph
       @method
   */
   this.shortestPath = shortestPathImpl;
   /**
       Test if the diagraph is DAG (Directed Acyclic Graph)
       @method
   */
   this.isDAG = isDAGImpl;
   /**
       Compute TC matrix if unweighed digraph, and distance matrix if weighted
       @method
   */
   this.warshallFloyd = warshallFloydImpl;
   /**
       Compute DFS-Based TC matrix
       @method
   */
   this.DfsTC = dfsTCImpl;
   /**
       Compute MST for undirected graph
       @method
   */
   //this.PrimMST = PrimImpl;
   this.PrimMST = primImpl2;
}

	// --------------------
	// more student methods next 




// -----------------------------------------------------------------------
// functions used by methods of Graph and ancillary objects

function addEdgeImpl(u_i, v_i)
{
    // fetch edge vertices using their id, where u: source vertex, v: target vertex
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

function addEdgeImpl2(u_i, v_i, weight)
{
    // fetch vertices using their id, where u: edge source vertex, v: target vertex
    var u = this.verts[u_i];
    var v = this.verts[v_i];

    // insert (u,v), i.e., insert v in adjacency list of u
    // (first create edge object using v_i as target, then pass object)
    var e = !(weight === undefined) ? new Edge(v_i, weight) : new Edge(v_i);
    u.adjacent.insert(e);

    // insert (v,u) if undirected graph (repeat above but reverse vertex order)
    if (!this.digraph)
    {
        e = !(weight === undefined) ? new Edge(u_i, weight) : new Edge(u_i);
        v.adjacent.insert(e);
    }

}

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
 
 
     // input vertex pairs from edge list input array
     // remember to pass vertex ids to add_edge()
 
     if (e.length > 0 && !(e[0].w === undefined))
         this.weighted = true;
     for (i = 0; i < this.ne; i++)
     {
         if (this.weighted)
             this.add_edge(e[i].u, e[i].v, e[i].w);
         else
             this.add_edge(e[i].u, e[i].v);
     }
 
     // double edge count if graph undirected
 
     if (!this.digraph)
     {
         this.ne = this.ne * 2;
     }


}

function printGraphImpl()
{
    document.write("<p>GRAPH {", this.label, "} ", this.weighted ? "WEIGHTED, " : "", this.digraph ? "" : "UN", "DIRECTED - ", this.nv, " VERTICES, ", this.ne, " EDGES:</p>");


    // report connectivity status if available

    document.write(this.componentInfo());


    // list vertices
    this.list_verts();
}

function dfsImpl(v_i)
{

    // get landing vert by id then process
    var v = this.verts[v_i];
    v.visit = true;
    len = this.dfs_push.length;
    this.dfs_push[len] = v_i;

    // recursively traverse unvisited adjacent vertices
    var w = v.adjacentByID();
    for (var j = 0; j < w.length; j++)
    {
        if (!this.verts[w[j]].visit)
        {
            this.dfs(w[j]);
        }
    }

}

// --------------------
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
        this.bfs_out[this.bfs_out.length] = u_i;


        // queue all unvisited vertices adjacent to u
        var w = u.adjacentByID();
        for (var i = 0; i < w.length; i++)
        {
            if (!this.verts[w[i]].visit) // mark and queue unvisited
            {
                this.verts[w[i]].visit = true;
                q.enqueue(w[i]);
            }
        }

    }

}
function topoSearchImpl(fun)
{
    nc = 0;

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
            fun == this.dfs ? this.dfs(i) : this.bfs(i);

            nc++;
        }
    }

    return nc;

}

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

// --------------------
function makeAdjMatrixImpl()
{
    // initially create row elements and zero the adjacency matrix
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


// Generate adjacency matrix representation of graph 
// obsolete, use makeAdjMatrixImpl3 instead. better encapsulation can be applied.

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
        var enodes = this.verts[i].incidentEdges();
        for (var j = 0; j < enodes.length; j++)
        {
            var edge_node = enodes[j];
            this.adjMatrix[i][edge_node.adjVert_i] = this.weighted ? edge_node.edgeWeight : 1;
        }
    }
}


function adjacentByIDImpl()
{
    var traversal = this.adjacent.traverse();

    var traversal_array = [];

    for (var i = 0; i < traversal.length; i++)
        traversal_array[i] = traversal[i].target_v;


    return traversal_array;

}

// Report connectivity information if available
function componentInfoImpl()
{
    return this.nc == 0 ? 'no connectivity info<br><br>': this.nc == 1 ? 'CONNECTED' : "DISCONNECTED: ".concat(this.nc, '<br>');
}


function insertAdjacentImpl(v_i, weight)
{
    this.adjacent.insert(new Edge(v_i, weight));
}

// 
function isConnectedImpl()
{
    this.nc = this.topoSearchImp(this.dfs); //Do topological search (DFS, or BFS) for the graph to set the number of connected components if not set
    return this.nc == 1; //The graph is connected if it has only one connected component
}

function vertexInfoImpl()
{
    var info = "".concat(" {", this.label, "} - VISIT: ", this.visit, " - ADJACENCY: ", this.adjacentByID(), "<br>");

    return info;
}

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


function makeGraphImpl(n, m, w) // no need for now
{

}
// -----------------------------------------------------------------------
// begin student code section
// -----------------------------------------------------------------------


// transitive closure package (CP8)

/**
	Check if there is a path between two vertices using their IDs
	@author Ghidaa Bajba
	@implements Graph#hasPath
	@param {integer} u_i Source vertex id
	@param {integer} v_i target vertex id
	@returns {boolean} True if there is path between u_i v_i
*/
function hasPathImpl(u_i, v_i)
{
    return this.warshallTC[u_i][v_i] == 1 ? true : false;
}

/**
	Return the shortest path between two vertices using their IDs
	@author Ghidaa Bajba
	@implements Graph#shortestPath
	@param {number} u_i Source vertex id
	@param {number} v_i target vertex id
	@returns {integer} The shortest path between u_i v_i
*/
function shortestPathImpl(u_i, v_i)
{
    return this.floydD[u_i][v_i];
}

/**
	Check if the given graph is Directed Acyclic Graph
	@author Ghidaa Bajba
	@implements Graph#isDAG
	@returns {boolean} True if diagraph is DAG
*/
function isDAGImpl()
{
    for (var i = 0, j = 0; i < this.warshallTC.length && j < this.warshallTC.length; i++, j++)
        if (this.hasPath(i, j))
            return false;
    return true;
}

/**
	Finds the shortest distance from one vertex to all the other vertices
	@author Ghidaa Bajba
	@implements Graph#warshallFloyd
*/
function warshallFloydImpl()
{
    // implement the ADJACENCY matrix 
    this.makeAdjMatrix();

    //Fill  warshallTC[] and distance matrices (floydD[]) by adjacent matrix
    for (var k = 0; k < this.adjMatrix.length; k++)
    {
        //Copy row by row
        this.warshallTC[k] = this.adjMatrix[k].slice();
        this.floydD[k] = this.adjMatrix[k].slice();
        for (var x = 0; x < this.nv; x++)
        {
            if (this.adjMatrix[k][x] == 0 && k != x)
            {
                this.floydD[k][x] = Infinity;
            }
        }
    }

    // warshall-Floyed algorithm
    for (var k = 0; k < this.floydD.length; k++)
    {
        for (var i = 0; i < this.floydD.length; i++)
        {
            for (var j = 0; j < this.floydD.length; j++)
            {
                this.floydD[i][j] = Math.min(this.floydD[i][j], (this.floydD[i][k] + this.floydD[k][j]));
                this.warshallTC[i][j] = this.warshallTC[i][j] || (this.warshallTC[i][k] && this.warshallTC[k][j]) ? 1 : 0;
            }
        }
    }

    //change the value from Infinity to 0 (because there is no distance = Infinity)
    for (var i = 0; i < this.floydD.length; i++)
        for (var j = 0; j < this.floydD.length; j++)
            if (this.floydD[i][j] == Infinity)
                this.floydD[i][j] = 0;

}

/**
	Finds the shortest distance from one vertex to all the other vertices (with better efficiency)
	@author Ghidaa Bajba
	@implements Graph#dfsTC
*/
function dfsTCImpl()
{
    // for each vertex
    for (var i = 0; i < this.nv; i++)
    {
        //process vertex v
        var v = this.verts[i];

        // mark all vertices unvisited
        for (var p = 0; p < this.nv; p++)
        {
            this.verts[p].visit = false;
        }

        // create and init the corresponding row 
        this.dfsTC[i] = [];
        for (var j = 0; j < this.nv; j++)
            this.dfsTC[i][j] = 0;

        //perform DFS search for each adjacent to the vertex v by its ID
        var w = v.adjacentByID();
        for (var n = 0; n < w.length; n++)
            this.dfs(w[n]); //for each adjacent vertex call dfs()

        //traverse the vertices to check which is visited
        for (var k = 0; k < this.nv; k++)
        {
            //if visited set 1 in the corresponding TC matrix
            if (this.verts[k].visit)
            {
                this.dfsTC[i][k] = 1;
            }
        }
    }
}

/**
	Greedy algorithm that finds a minimum spanning tree for a weighted undirected graph.
	@author Ghidaa Bajba
	@implements Graph#PrimMST
*/
function PrimImpl()
{
    // create set of vertices tree 
    var n;
    var verticesTree = [];

    // make all vertex unvisited
    for (var m = 0; m < this.nv; m++)
    {
        this.verts[m].visit = false;
    }
    // inatiat first index with first value on vert array 
    verticesTree[0] = this.verts[0];
    this.verts[0].visit = true;

    var min = Infinity; // any number would be less than infinty 
    for (var i = 0; i < this.nv; i++)
    {
        for (var j = 0; j < verticesTree.length; j++)
        {
            //get fringe vertices of the current vertext
            var incident_Edge = verticesTree[j].incidentEdges();
            //loop on every fringe vertex
            for (var k = 0; k < incident_Edge.length; k++)
            {
                //check every unvisited vertex to check if it can be visited in a shorter distance
                if ((!this.verts[incident_Edge[k].adjVert_i].visit) && (incident_Edge[k].edgeWeight < min))
                {
                    this.Prim_Edge[i] =
                        (
                        {
                            v: verticesTree[j],
                            u: this.verts[incident_Edge[k].adjVert_i],
                            w: incident_Edge[k].edgeWeight
                        });
                    //set the new weight as the minimum
                    min = this.Prim_Edge[i].w;
                }
            }
        }
        n = this.Prim_Edge.length;
        verticesTree[verticesTree.length] = this.Prim_Edge[n - 1].u;

        // mark VerticesTree as visited
        this.Prim_Edge[n - 1].u.visit = true;

        min = Infinity;
    }
}

// -----------------------------------------------------------------------
// // published docs section (ref. specs, assignment page)
// use starter6-based P1M1 code as-is (fixes/improvements OK)
// no JSDOC comments in this section (docs already published)
// -----------------------------------------------------------------------

function list_verts()
{
	var i, v;  
	for (i=0; i < this.nv; i++)
	{
	   v = this.verts[i];
	   document.write( "VERTEX: ", i, " {", v.label, "} - VISIT: ", v.visit,
		  " - ADJACENCY: ", v.adjacentByID(), "<br>" );
	}
 
} 

/**
 * * prim is an Greedy algorithm that finds a minimum spanning tree for a weighted undirected graph. minimal spanning tree (MSP). 
 * It finds the connected acyclic subgraph (tree) that contains all vertices with smallest total weight.
 * 
   @author Ghidaa Bajba
   @implements Graph#PrimMST
*/
function primImpl2(){

  
    //chech if the graph is weighted, otherwise there is no minimal spanning tree. 
    
    if(this.weighted){
   
       var pq = new PQueue();  // create qeueue 
   
       var penultimate = []; // array contains penulimate
   
       var weight = []; //  array contains weight of edges
   
       this.VT=[]; //create an array of json objects (it's minimum spanning tree)
   
       var u;
   
   
   
       // mark all vertices unvisited
       // intitalize penulimate array by '-' that mean's no penulimate until now.
       // intitalize weight array by infinity.
   
       for (var j = 0; j < this.nv; j++)
   
       {
   
           this.vert[j].visit = false;
   
           penultimate[j] = "-";
   
           weight[j] = Infinity;
   
       }
   
   
   
       // edges array contains all the edges (v,u)
       var edges =this.vert[0].incidentEdges();
   
       for(var i=0;i<edges.length;i++){
   
           //insert dest. vertex and weight of edge in the queue.
           pq.insert(edges[i].adjVert_i, edges[i].edgeWeight); 
   
           // make a src. vertex penultimate of dest. vertex.
           penultimate[edges[i].adjVert_i]=0; 
   
           //add weight of edge to weight array.
           weight[edges[i].adjVert_i] = edges[i].edgeWeight;
   
       }
   
   
   
       //intitalize Vt by the first vertex and mark it as true.
       this.vert[0].visit = true;
   
       this.VT[0] = {
   
           //vertex
           tree: 0,
   
           //penultimate of this vertex
           parent: penultimate[0]
   
       } 
   
   
   
       //find a way to the rest n-1 vertices
   
       for(var i=1; i<this.nv; i++){
   
           //delete minimum-weight edge from queue
           u =pq.deleteMin();
   
           //mark a vertex as visit
           this.vert[u].visit = true;
   
          
           //add vertex that delete from qeueue into vt
           this.VT[i] = {
   
               tree: u,
   
               parent: penultimate[u]
   
           }  
   
   
   
          // edges array contains all the edges (v,u)
           edges=this.vert[u].incidentEdges();
   
           for(var j=0;j<edges.length;j++){
   
   
            //if vertex doesn't mark as visit which mean's don't add into VT & weight of edge less than weight of edge in qeueue 
            //then insert or update  dest. vertex and its weight ,
           //update penultimate and update weight of edge in weight array. 
               if(!this.vert[edges[j].adjVert_i].visit && edges[j].edgeWeight < weight[edges[j].adjVert_i]){
   
                   pq.insert(edges[j].adjVert_i,edges[j].edgeWeight);
   
                   penultimate[edges[j].adjVert_i]=u;
   
                   weight[edges[j].adjVert_i] = edges[j].edgeWeight;
   
               }
               
   
           }
   
       }
   
    }
   
   }
   
   
/**
 * Dijkstra is an algorithm for single-source-shortest-paths problem. 
 * this algorithm is applicable to undirected and directed graphs with nonnegative weigths.
 * this algorithm works for weighted connectd graph
 * It finds the shortest paths to a graph's vertices in order of their distanec from a given source.
 * It represnts row s in the distance matrix
 * the tree edges will be stored in the global variable VT
 * 
   @author Hadeel ALoufi
   @implements Graph#Dijkstra
   @param {integer} s the source vertex
*/

function DijkstraImpl(s)
{
//chech if the graph is weighted, otherwise the shortest path will not be computed
if(this.weighted){

// create a queue
var pq = new PQueue();

// create arrays
var distances_v = [];
var penultimate_v = [];

//initialize queue
for(var i=0; i<this.nv; i++)
{
    //distances[i] is the shortest path length from s to v
    distances_v[i] = Infinity;

    //penultimate[i] is the next-to-last vertex in path
    penultimate_v[i] = "-";

    //insert in queue
    pq.insert(i, distances_v[i]);

    //mark all vertices as unvisited. this meaning is not in the tree
    this.verts[i].visit = false;
}

//initialize the tree by the (s) source
distances_v[s] = 0;
pq.insert(s, distances_v[s]);

//get the next vertex to be added to the tree
for(var i=0; i<this.nv; i++)
{
    //pick a fringe vertex with min distance
    u = pq.deleteMin();

    //mark all vertices as visited.
    this.verts[u].visit = true;

    //the tree edges will be stored in VT
    this.VT[i] = 
    {
        //tree represents the vertex.
        tree: u,
        //parent represents the next-to-last vertex in path.
        parent: penultimate_v[u],
        //distance represents the total distance from the source (s) to the current vertex.
        distance: distances_v[u]
    }

    //get adjacents of the new added vertex
    var adjacents = this.verts[u].incidentEdges();

    //update fringe set after adding fringe vertex
    for(var k=0; k<adjacents.length; k++)
    {
        var id = adjacents[k].adjVert_i;
        var weight = adjacents[k].edgeWeight;

        //if there's a shorter path to a vertex, consider that path (update the queue, distances, and penultimate arrays) 
        if(distances_v[u] + weight < distances_v[id]&&!this.verts[id].visit){
                distances_v[id] = distances_v[u] + weight;
                penultimate_v[id] = u;
                pq.insert(id,distances_v[id]);
        } 
    }
}
}
}
function printTree(VT){
    for(i=0; i<VT.length-1; i++){
        document.write((VT[i].distance===undefined? '':VT[i].distance) ,'(',VT[i].parent,',',VT[i].tree,'),');
    }
    document.write((VT[VT.length-1].distance===undefined? '':VT[VT.length-1].distance) ,'(',VT[VT.length-1].parent,',',VT[VT.length-1].tree,').<br>');
}
