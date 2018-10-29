// pad a string with 0's on the left
function padLeft(str,len){
    return (Array(len).join("0") + str).slice(-len);
}

function base_p_expansion(n, p, padlength = 0){
    let result = [];
    while(n!==0){
        result.push(n % p);
        n = (n/p)|0;
    }
    while(result.length < padlength ){
        result.push(0);
    }
    return result;
}

// Integral binomial coefficient.
function binomial(n, k) {
    if ((typeof n !== 'number') || (typeof k !== 'number'))
        return false;
    let coeff = 1;
    for (let x = n-k+1; x <= n; x++) coeff *= x;
    for (let x = 1; x <= k; x++) coeff /= x;
    return coeff;
}


// Mod 2 multinomial coefficient
function multinomial_mod_2(list){
    let old_sum = list[0];
    for(let i = 1; i < list.length; i++){
        for(let j = 1; j <= Math.min(old_sum, list[i]); j = j<<1) {
            if (((j & old_sum ) === j) && ((j & list[i]) !== 0) ) {
                return 0;
            }
        }
        old_sum += list[i];
    }
    return 1;
}

function binomial_mod2(n,k){
    if(n < k){
        return 0
    } else {
        return +!!(((n-k) & k) === 0);
    }
}


function multinomial_odd(list,p){
    let n = list.reduce((a,b)=>a+b);
    let answer = 1;
    let n_expansion = base_p_expansion(n, p);
    let list_expansion = list.map(x => base_p_expansion(x,p,n_expansion.length));
    let index = 0;
    for(let index = 0; index < n_expansion.length; index++) {
        let multi = 1;
        let partial_sum = 0;
        for(let exp of list_expansion) {
            if(index < exp.length){
                partial_sum = partial_sum + (+exp[index]);
                multi = (multi * binomial(partial_sum, +exp[index]));
            }
        }
        answer = (answer * multi) % p;
        if(answer == 0){
            return 0;
        }
    }
    return answer
}

function binomial_modp(n,k,p){
    if(n < k){
        return 0;
    }
    return multinomial_odd([n-k, k], p);
}

function multinomial(list, p = 2){
    if(p===2){
        return multinomial_mod_2(list);
    } else {
        return multinomial_odd(list, p);
    }
}

function binomial_gen(n,k,p = 2){
    if(n<k || k<0){
        return 0;
    }
    if(p===2){
        return binomial_mod2(n,k);
    } else {
        return binomial_modp(n,k,p);
    }
}

exports.multinomial = multinomial;
exports.binomial = binomial_gen;