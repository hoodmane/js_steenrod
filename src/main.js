let stringify = require("json-stringify-pretty-compact")
let preprocess = require("./preprocessor");

let SerreCartanBasis = require("./SerreCartanBasis");
let MilnorBasis = require("./MilnorBasis");
require("./ChangeBasis");

context = {};
context.basis = SerreCartanBasis;
context.p = 2;

function passToBasis(f){
    context[f + "_call"] = function(...args){
        return context.basis[f](args,context.p);
    }
}

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

}

context.MINUS = function(a){
    if(a.scale){
        return a.scale(-1);
    } else {
        return -a;
    }
};



["Sq", "b", "P", "Q", "pst"].forEach(passToBasis);

function run(code){
    let result = eval(preprocess(code));
    result = stringify(result);
    if(result){
        result = result.replace(/"1s@/g,"").replace(/@s1"/g, "");
    }
    return result;
}

//console.log(run("Sq2*Sq1"));

if(typeof window !== 'undefined'){
    window.context = context;
    window.run = run;
}