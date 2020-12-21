// graph from Exercise 10.2 2a

var _v2 = [
    {label:"1"}, //index 0
    {label:"2"}, //index 1
    {label:"3"}, //index 2
    {label:"4"}, //index 3
    {label:"5"}, //index 4
    {label:"6"}  //index 5
]; 

var _e2 = [
    // 1(0) --> 2(1)
    {u:0, v:1, w:2},
    
    // 1(0) --> 3(2)
    {u:0, v:2, w:7},
    
    // 2(1) --> 4(3)
    {u:1, v:3, w:3},

    // 2(1) --> 5(4)
    {u:1, v:4, w:4},

    // 3(2) --> 4(3)
    {u:2, v:3, w:4},

    // 4(3) --> 6(5)
    {u:3, v:5, w:1},

    // 5(4) --> 6(5)
    {u:4, v:5, w:5},
]