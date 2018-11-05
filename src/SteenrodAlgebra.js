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

class BasisWrapper {
    constructor(alg, v){
        this.algebra = alg;
        if(v === undefined){
            return;
        }
        if(v.constructor === MilnorBasis){
            this.milnor = v;
        } else if(v.constructor === SerreCartanBasis){
            this.serre_cartan = v;
        }
    }

    get(basis){
        if(basis === SerreCartanBasis){
            return this.getSerreCartan();
        }
        if(basis === MilnorBasis) {
            return this.getMilnor()
        }
        throw new Error("Unknown basis");
    }

    getMilnor(){
        if(!this.milnor){
            this.milnor = this.serre_cartan.toMilnor();
        }
        return this.milnor;
    }

    getSerreCartan(){
        if(!this.serre_cartan){
            this.serre_cartan = this.milnor.toSerreCartan();
        }
        return this.serre_cartan;
    }

    degree(){
        if(this.serre_cartan){
            return this.serre_cartan.degree();
        } else {
            return this.milnor.degree();
        }
    }

    excess(){
        if(this.serre_cartan){
            return this.serre_cartan.excess();
        } else {
            return this.milnor.excess();
        }
    }

    antipode(){
        let sc = this.getSerreCartan();
        throw new Error("Not implemented");
    }

    add(other_vector){
        if(this.algebra.p !== otherVector.algebra.p){
            throw new Error(`Vectors ${this} and ${other_vector} belong to different algebras and can't be multiplied.`);
        }
        let result = new BasisWrapper(this.algebra);
        if(this.serre_cartan && other_vector.serre_cartan){
            result.serre_cartan = this.serre_cartan.add(other_vector.serre_cartan)
        }
        if(this.milnor && other_vector.milnor){
            result.milnor = this.milnor.add(other_vector.milnor);
        }
        if(!result.serre_cartan && !result.milnor){
            result.serre_cartan = this.getSerreCartan().add(other_vector.getSerreCartan());
        }
        return result;
    }

    mult(other_vector){
        if(this.algebra.p !== other_vector.algebra.p){
            throw new Error(`Vectors ${this} and ${other_vector} belong to different algebras and can't be multiplied.`);
        }
        let result = new BasisWrapper(this.algebra);
        if(this.serre_cartan && other_vector.serre_cartan){
            console.log("sc: ", this.serre_cartan);
            console.log("sc: ",other_vector.serre_cartan);
            result.serre_cartan = this.serre_cartan.mult(other_vector.serre_cartan)
        }
        if(this.milnor && other_vector.milnor){
            result.milnor = this.milnor.mult(other_vector.milnor);
        }
        if(!result.serre_cartan && !result.milnor){
            result.serre_cartan = this.getSerreCartan().mult(other_vector.getSerreCartan());
        }
        return result;
    }

    equal(otherVector){
        if(this.algebra.p !== otherVector.algebra.p){
            return false;
        }
        if(this.serre_cartan){
            return this.serre_cartan.equal(otherVector.getSerreCartan());
        } else {
            return this.milnor.equal(otherVector.getMilnor());
        }
    }

    toString(){
        return this.get(this.algebra.basis).toString();
    }

    toJSON(){
        return "1s@" + this.toString() + "@s1";
    }
}

BasisWrapper.basis = SerreCartanBasis;

class SteenrodAlgebra {
    constructor(prime, basis){
        this.p = prime;
        this.basis = basis;
    }

    Sq(...n){
        if(n.length > 1){
            return new BasisWrapper(this, MilnorBasis.Sq(n, this.p));
        } else {
            return new BasisWrapper(this, this.basis.Sq(n, this.p));
        }
    }

    b(){
        return new BasisWrapper(this, this.basis.b());
    }

    P(...n){
        if(n.length > 1){
            return new BasisWrapper(this, MilnorBasis.P(n, this.p));
        } else {
            return new BasisWrapper(this, this.basis.P(n, this.p));
        }
    }

    bP(...n){
        if(n.length > 1){
            return new BasisWrapper(this, MilnorBasis.bP(n, this.p));
        } else {
            return new BasisWrapper(this, this.basis.bP(n, this.p));
        }
    }

    Q(...n){
        return new BasisWrapper(this, MilnorBasis.Q(n, this.p));
    }

    pst(...st){
        console.log(st);
        return new BasisWrapper(this, MilnorBasis.pst(st, this.p));
    }
}



module.exports = SteenrodAlgebra;




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