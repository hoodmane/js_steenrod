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
        return this.get(this.algebra.getBasis()).toString();
    }

    toJSON(){
        return "1s@" + this.toString() + "@s1";
    }
}

BasisWrapper.basis = SerreCartanBasis;

class SteenrodAlgebra {
    constructor(prime, basis){
        this.p = prime;
        this._generic = this.p !== 2;
        this._basis = basis;
    }

    prime(){
        return this.p;
    }

    one(){
        return this.basis(0)[0];
    }

    setBasis(b){
        this._basis = b;
    }

    Sq(...n){
        if(n.length > 1){
            return new BasisWrapper(this, MilnorBasis.Sq(n, this.p));
        } else {
            return new BasisWrapper(this, this._basis.Sq(n, this.p));
        }
    }

    b(){
        return new BasisWrapper(this, this._basis.b());
    }

    P(...n){
        if(n.length > 1){
            return new BasisWrapper(this, MilnorBasis.P(n, this.p));
        } else {
            return new BasisWrapper(this, this._basis.P(n, this.p));
        }
    }

    bP(...n){
        if(n.length > 1){
            return new BasisWrapper(this, MilnorBasis.bP(n, this.p));
        } else {
            return new BasisWrapper(this, this._basis.bP(n, this.p));
        }
    }

    Q(...n){
        return new BasisWrapper(this, MilnorBasis.Q(n, this.p));
    }

    pst(...st){
        return new BasisWrapper(this, MilnorBasis.pst(st, this.p));
    }

    getBasis(){
        return this._basis;
    }

    basis(n){
        let result = [];
        for(let basis_elt of this._basis.basis(n, this.p, this.generic)){
            let vec = new this._basis(this.p);
            vec.set(basis_elt, 1);
            result.push(new BasisWrapper(this, vec));
        }
        return result;
    }

    antipode_on_basis(t){
        /*
        The antipode of a basis element of this algebra

        INPUT:

        - ``t`` -- tuple, the index of a basis element of self

        OUTPUT: the antipode of the corresponding basis element,
        as an element of self.

        ALGORITHM: according to a result of Milnor's, the antipode of
        `\text{Sq}(n)` is the sum of all of the Milnor basis elements
        in dimension `n`. So: convert the element to the Serre-Cartan
        basis, thus writing it as a sum of products of elements
        `\text{Sq}(n)`, and use Milnor's formula for the antipode of
        `\text{Sq}(n)`, together with the fact that the antipode is an
        antihomomorphism: if we call the antipode `c`, then `c(ab) =
        c(b) c(a)`.

        At odd primes, a similar method is used: the antipode of
        `P(n)` is the sum of the Milnor P basis elements in dimension
        `n*2(p-1)`, multiplied by `(-1)^n`, and the antipode of `\beta
        = Q_0` is `-Q_0`. So convert to the Serre-Cartan basis, as in
        the `p=2` case.

        EXAMPLES::

            sage: A = SteenrodAlgebra()
            sage: A.antipode_on_basis((4,))
            Sq(1,1) + Sq(4)
            sage: A.Sq(4).antipode()
            Sq(1,1) + Sq(4)
            sage: Adem = SteenrodAlgebra(basis='adem')
            sage: Adem.Sq(4).antipode()
            Sq^3 Sq^1 + Sq^4
            sage: SteenrodAlgebra(basis='pst').Sq(3).antipode()
            P^0_1 P^1_1 + P^0_2
            sage: a = SteenrodAlgebra(basis='wall_long').Sq(10)
            sage: a.antipode()
            Sq^1 Sq^2 Sq^4 Sq^1 Sq^2 + Sq^2 Sq^4 Sq^1 Sq^2 Sq^1 + Sq^8 Sq^2
            sage: a.antipode().antipode() == a
            True

            sage: SteenrodAlgebra(p=3).P(6).antipode()
            P(2,1) + P(6)
            sage: SteenrodAlgebra(p=3).P(6).antipode().antipode()
            P(6)

        TESTS::

            sage: Milnor = SteenrodAlgebra()
            sage: all([x.antipode().antipode() == x for x in Milnor.basis(11)]) // long time
            True
            sage: A5 = SteenrodAlgebra(p=5, basis='adem')
            sage: all([x.antipode().antipode() == x for x in A5.basis(25)])
            True
            sage: H = SteenrodAlgebra(profile=[2,2,1])
            sage: H.Sq(1,2).antipode() in H
            True
        */
        let p = this.prime();
        // if(this.getBasis() !== SerreCartanBasis){
        //     return this._change_basis_on_basis(t, 'serre-cartan').antipode();
        // }
        let antipode = this.one();
        if(!this._generic) {
            for (let n of t) {
                antipode = self(sum(SteenrodAlgebra().basis(n))) * antipode
            }
            return antipode
        }
        let index = 0;
        for(let n of t){
            if( index % 2 === 0 ){
                if(n !== 0){
                    antipode = -self.Q(0) * antipode
                }
            } else {
                B = SteenrodAlgebra(p=p,generic=self._generic).basis(n * 2 * (p-1))
                s = self(0)
                for(let b of B){
                    if(len(b.leading_support()[0]) == 0){
                        s += self(b)
                    }
                }
                antipode = (-1)**n * s * antipode;
            }
        }
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