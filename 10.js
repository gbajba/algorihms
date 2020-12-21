// CPCS 324 Algorithms & Data Structures 2
// Outline - Priority queue data structure
// 2019, Dr. Muhammad Al-Hashimi


// -----------------------------------------------------------------------
// Basic design decisions and implementation planning (objects & interfaces)

// initial requirements: to quickly support Dijkstra and second Prim's algorithms, 
// implement minimum priority functionality

// design decisions:
// Reuse the 324 linked list implementation
// how will the PQ ops be implemented?
// <student fill here>

// we designed a priority queue that based on linked list package.
// A priority queue needs a data-item and a priority key. Threrefore,
// we should pass these tow things as an object to the linked list which
// is designed to store one item.
// We  decided to use the choice number 4 of the four choices listed in
// (First PQ Design Notes) page.
// In this choice the PQ methods of PQ operations can directly process the linked list.
// Therefore, from the PQ package we can directly manipulate linked list internals. 

// code plan: start to write your API specifications (JSDOC comments) and think about 
// design consequences (design impact)

// Impact analysis:
// <student fill here>

// As a result of the choosen design choice, we will not maintain a 
// proper encapsulation and the choice violates the principles of
// good OOP (object-oriented programming). We used three mthods from
// linked list package. Three methods are : isEmpty(), deleteFirst() 
// and insert().
// Time complexity for each method : 
// - isEmpty() : O(1)
// - deleteMin : O(n)
// - insert() : O(1) for insert case, or O(n) for update case 




// -----------------------------------------------------------------------
// Priority queue object constructor (document using JSDOC comments)

/**
   Create a priority queue

   @author Asma Megabel

   @constructor
*/
function PQueue()
{
    /** Head pointer of a queue  */
    this.pq = new List();

    // specify (design) methods
    /**
      Return true if queue is empty
      @method
    */
    this.isEmpty = isEmptyImpl;

    /**
      Remove and return item with minimum priority
      @method
     */
    this.deleteMin = deleteMinImpl;

    /**
      Insert or update an item with priority
      @method
     */
    this.insert = insertImpl;

}

// -----------------------------------------------------------------------
// Priority queue node constructor (document using JSDOC comments)

/**
   Create a new PQ node with two parameters item and its key priority

   @author Asma Megabel

   @constructor
   @param {integer} item Data item value (vertex id)
   @param {integer} key Priority key value 
*/

function PQNode(item, key)
{
    /**
      The value of data item 
    */
    this.item = item;

    /**
      The value of priority 
    */
    this.prior = key;

    // specify (design) methods

}

// -----------------------------------------------------------------------
// functions used by PQueue() object methods
// specify interface information (JSDOC comments)
// function names should not clash with linklist.js and queue.js
// ....

/**
   Return true if the queue is empty otherwise return false.
   We used isEmpty() method from linked list package.
   
   @author Asma Megabel

   @implements PQueue#isEmpty
   @returns {boolean} True if PQ is empty 
*/


function isEmptyImpl()
{

}

/**
   Remove the item with the minimum priority from 
   the priority queue and return the data item 
   of removed item.

   Traverse the priority queue for searching for the item with the minimum
   value of priority key. If that item is the first one in priority queue then 
   delete it using deleteFirst() method from linked list package otherwise change
   pointers to delete the item.

   @author Asma Megabel

   @implements PQueue#deleteMin
   @return {integer} The data item of deleted item with the minimum priority
*/


function deleteMinImpl()
{

}

/**
   Insert a new item into the priority queue
   or update an item with new priority if the new key 
   smaller than the current key.

   Insert the new item if the priority queue is empty otherwise traverse the
   priority queue to check if the item exist in the queue or not. If the item is
   exist then it is the update case otherwise it is the insert case. We used 
   insert() method from linked list package in insert case.

   @author Asma Megabel

   @implements PQueue#insert
   @param {ineger} item Data item to be inserted in the PQ
   @param {integer} key Priority key value  
*/

function insertImpl(item, key)
{

}
