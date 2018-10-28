let StringifyingMap = require("./StringifyingMap");

class Vector extends StringifyingMap {
    constructor(p){
        super(JSON.stringify);
        this.p = p;
    }

    addTerm(monomial, coeff){
        if(this.has(monomial)){
            console.log(this.get(monomial));
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

    toString(){
        return this._m;
    }
}

module.exports = Vector;