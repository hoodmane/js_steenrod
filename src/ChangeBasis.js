let SerreCartanBasis = require("./SerreCartanBasis");
let MilnorBasis = require("./MilnorBasis");

SerreCartanBasis.basis_to_milnor = function(b, p){
    let result = MilnorBasis.unit(p);
    if(p===2){
        for(let j of b){
            result = result.mult(MilnorBasis.Sq(j,p));
        }
    } else {
        let i = 0;
        for(let j of b){
            i++;
            if(i%2 === 1){
                if(j!==0){
                    result = result.mult(MilnorBasis.Q(0,p));
                }
            } else {
                result = result.mult(MilnorBasis.P(j,p))
            }
        }
    }
    return result;
};

SerreCartanBasis.prototype.toMilnor = function(){
    return this.map_on_basis((x) => SerreCartanBasis.basis_to_milnor(x,this.p),MilnorBasis);
};

let pad_array = function(arr, len, fill) {
    return arr.concat(Array(len).fill(fill)).slice(0,len);
};

MilnorBasis.basis_to_serrecartan = function(b, p) {
    let result;
    if(p === 2) {
        // See Monks paper page 8
        let t = Array(b.length);
        t[b.length - 1] = b[b.length - 1];
        for (let i = b.length - 2; i >= 0; i--) {
            t[i] = b[i] + 2 * t[i + 1];
        }
        let x = SerreCartanBasis.basis_to_milnor(t, p);
        x.delete(b);
        result = x.map_on_basis((m) => MilnorBasis.basis_to_serrecartan(m, p),SerreCartanBasis);
        result.set(t,1);
    } else {
        let e = b[0];
        let s = b[1];
        let len = Math.max(s.length, ...e);
        s = pad_array(s, len, 0);
        let t = Array(2*len + 1);
        t.fill(0);
        for(let i of e){
            t[2*i] = 1;
        }
        t[t.length - 2] = s[s.length - 1] + t[t.length - 1];
        let idx = t.length - 2;
        for (let i = s.length - 2; i >= 0; i--) {
            idx -= 2;
            t[idx] = t[idx + 1] + s[i] + p * t[idx + 2];
        }
        let x = SerreCartanBasis.basis_to_milnor(t, p);
        x.delete(b);
        x.scale(-1);
        result = x.map_on_basis((m) => MilnorBasis.basis_to_serrecartan(m, p),SerreCartanBasis);
        result.set(t,1);
    }
    let result_copy = new SerreCartanBasis(p);
    for(let [k,v] of result){
        result_copy.set(k,v);
    }
    return result_copy;
};

MilnorBasis.prototype.toSerreCartan = function(){
    return this.map_on_basis((x) => MilnorBasis.basis_to_serrecartan(x,this.p),SerreCartanBasis);
};

// console.log("a:");
// console.log(SerreCartanBasis.basis_to_milnor([0,1,1],3));
// console.log(SerreCartanBasis.basis_to_milnor([2,1],2));
// console.log(SerreCartanBasis.basis_to_milnor([0,1,1],3));
// console.log(MilnorBasis.basis_to_serrecartan([0,1],2));
// console.log(MilnorBasis.basis_to_serrecartan([1,1],2));
// console.log(MilnorBasis.basis_to_serrecartan([0,2],2));
// console.log(MilnorBasis.basis_to_serrecartan([1,2],2));
// console.log("----");
//
// console.log(MilnorBasis.basis_to_serrecartan([[],[0,1]],3));
// console.log(MilnorBasis.basis_to_serrecartan([[],[1,1]],3));
// console.log(MilnorBasis.basis_to_serrecartan([[],[0,2]],3));
// console.log(MilnorBasis.basis_to_serrecartan([[],[1,2]],3));
// console.log(MilnorBasis.basis_to_serrecartan([[],[0,3]],3));
// console.log(MilnorBasis.basis_to_serrecartan([[],[1,3]],3));
//console.log(MilnorBasis.basis_to_serrecartan([[1],[0,1]],3));


// console.log(SerreCartanBasis.basis_to_milnor([0,3,1],3));
// console.log(SerreCartanBasis.basis_to_milnor([0,3,0,1,0],3));
// console.log(SerreCartanBasis.basis_to_milnor([2,1],2));
// console.log(SerreCartanBasis.basis_to_milnor([3,1],2));
// console.log(SerreCartanBasis.basis_to_milnor([4,1],2));
// console.log(SerreCartanBasis.basis_to_milnor([4,2],2));