// graph from {Figure 10.4 (Levitin, 3rd edition)}

var _v = [
    { label: "1" },	//Index = 0
    { label: "2" },	//Index = 1
    { label: "3" },	//Index = 2
    { label: "4" },	//Index = 3
    { label: "5" },	//Index = 4
    { label: "6" },	//Index = 5
    ];
    
    //There are 6 verticese and 7 edges
    var _e = [
    // 1(0) --> 2(1)
    { u: 0, v: 1, w: 2 },

    // 1(0) --> 4(3)
    { u: 0, v: 3, w: 3 },

    // 2(1) --> 3(2)
    { u: 1, v: 2, w: 5 },

    // 2(1) --> 5(4)
    { u: 1, v: 4, w: 3 },

    // 3(2) --> 6(5)
    { u: 2, v: 5, w: 2 },
    
    // 4(3) --> 3(2)
    { u: 3, v: 2, w: 1 },
    
    // 5(4) --> 6(5)
    { u: 4, v: 5, w: 4 },
    
    ];