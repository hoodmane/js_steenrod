let StringifyingMap = require("./StringifyingMap");

class Vector extends StringifyingMap {
    constructor(p){
        super(JSON.stringify);
        this.p = p;
    }

    copy(){
        let result = new this.constructor(this.p);
        for(const [mono,coeff] of this){
            result.set(mono,coeff);
        }
        return result;
    }

    addTerm(monomial, coeff){
        if(this.has(monomial)){
            coeff += this.get(monomial);
        }
        coeff = (coeff % this.p + this.p) % this.p;
        if(coeff === 0){
            this.delete(monomial);
        } else {
            this.set(monomial, coeff);
        }
        return this;
    }

    addVector(vect, coeff = 1){
        for(const [key, value] of vect){
            this.addTerm(key, coeff*value);
        }
        return this;
    }

    scale(s){
        for(const [key, value] of this){
            this.set(key, s*value);
        }
        return this;
    }

    add(other_vector){
        let result = new this.constructor(this.p);
        for(const [mono,coeff] of this){
            result.set(mono,coeff);
        }
        for(const [mono,coeff] of other_vector){
            result.addTerm(mono,coeff);
        }
        return result;
    }

    mult(term_multiply, other_vector){
        // This method of getting an object of the same time as the subclass is a bit limited...
        let result = new this.constructor(this.p);
        for(const [mono1,coeff1] of this){
            for(const [mono2,coeff2] of other_vector){
                let prod = term_multiply(mono1,mono2);
                result.addVector(prod,coeff1*coeff2);
            }
        }
        return result;
    }

    equal(other_vector){
        let checkedKeys = new StringifyingMap();
        for(let [k, v] of this){
            if(this.get(k) !== other_vector.get(k)){
                return false;
            }
            checkedKeys.set(k,1);
        }
        for(let [k, v] of other_vector){
            if(!checkedKeys.has(k)){
                return false;
            }
        }
        return true;
    }

    map_on_basis(f, output_type){
        if(output_type === undefined){
            output_type = this.constructor;
        }
        let result = new output_type(this.p);
        for(const [mono1,coeff1] of this){
            result.addVector(f(mono1),coeff1);
        }
        return result;
    }
}

module.exports = Vector;