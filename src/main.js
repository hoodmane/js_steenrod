let stringify = require("json-stringify-pretty-compact")
let preprocess = require("./preprocessor");

let SerreCartanBasis = require("./SerreCartanBasis");
let MilnorBasis = require("./MilnorBasis");
let SteenrodAlgebra = require("./SteenrodAlgebra");

let algebra = new SteenrodAlgebra(2, SerreCartanBasis);
context = {};
context.algebra = algebra;

function setPrime(p){
    context.algebra = new SteenrodAlgebra(p, context.algebra.basis);
}

bases = {"Milnor" : MilnorBasis, "SerreCartan" : SerreCartanBasis};

function setBasis(b){
    if(bases[b]){
        b = bases[b];
    }
    context.algebra.basis = b;
}


setBasis(SerreCartanBasis);

context.MULT = function(a,b){
    if(a.mult){
        return a.mult(b);
    } else if(b.mult){
        return b.mult(a);
    } else {
        return a*b;
    }
};

context.ADD = function(a,b){
    if(a.add){
        return a.add(b);
    } else if(b.add){
        return b.add(a);
    } else {
        return a + b;
    }
};

context.SUB = function(a,b){
    if(a.add && b.scale){
        return a.add(b.scale(-1));
    } else {
        return a - b;
    }

};

context.MINUS = function(a){
    if(a.scale){
        return a.scale(-1);
    } else {
        return -a;
    }
};

context.EQUAL = function(x, y){
    if(x.equal){
        return x.equal(y);
    }
    if(y.equal){
        return false;
    }
    return x == y;
};


function passToAlgebra(f){
    context[f] = function(...args){
        return context.algebra[f](...args);
    }
}

["Sq", "b", "P","bP", "Q", "pst"].forEach(passToAlgebra);


function run(code){
    let result = eval(preprocess(code));
    result = stringify(result);
    if(result){
        result = result.replace(/"1s@/g,"").replace(/@s1"/g, "");
    }
    return result;
}

// console.log("sq2: ", context.Sq(2).mult(context.Sq(1)).getSerreCartan());
//
console.log(run("Sq2*Sq1"));
// console.log(run("Sq1==Sq1"));
// console.log(run("b"));
// console.log(SerreCartanBasis.b(2));


//console.log(context.Q(1,1));
//console.log(run("range(4).map(Sq)"));

function range(start, stop, step = 1){
    if(arguments.length === 1){
        start = 1;
        stop = arguments[0];
        step = 1;
    }
    return Array(Math.ceil((stop - start + step)/step)).fill(start).map((x, y) => x + y * step);
}


if(typeof window !== 'undefined'){
    window.context = context;
    window.setPrime = setPrime;
    window.setBasis = setBasis;
    window.SteenrodAlgebra = SteenrodAlgebra;
    window.run = run;
    window.range = range;
}