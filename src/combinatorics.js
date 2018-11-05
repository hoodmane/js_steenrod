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


function restricted_partitions(n, l, no_repeats=false){
    /*
            restricted_partitions(10, [6,4,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
            restricted_partitions(10, [6,4,4,4,2,2,2,2,2,2])
            [[6, 4], [6, 2, 2], [4, 4, 2], [4, 2, 2, 2], [2, 2, 2, 2, 2]]
    */
    if(n < 0){
        return [];
    } else if(n === 0){
        return [[]];
    } else {
        let results = [];
        let index;
        if(no_repeats){
            index = 1;
        } else {
            index = 0;
        }
        let old_i = 0;
        for(let i of l){
            if(old_i !== i){
                for(let sigma of restricted_partitions(n-i, l.slice(index), no_repeats)){
                    results.push([i, ...sigma]);
                }
            }
            index ++;
            old_i = i;
        }
        return results;
    }
}

function xi_degrees(n,p=2, reverse=true){
    /*
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17)
            [15, 7, 3, 1]
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17, reverse=False)
            [1, 3, 7, 15]
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(17,p=3)
            [13, 4, 1]
            sage: sage.algebras.steenrod.steenrod_algebra_bases.xi_degrees(400,p=17)
            [307, 18, 1]
    */
    if(n <= 0){
        return []
    }
    let N = n*(p-1) + 1;
    let xi_max = 1;
    while( N > 0){
        N = Math.floor(N/p);
        xi_max ++;
    }
    let l = [];
    for(let d = 1; d < xi_max; d++){
        l.push(((Math.pow(p,d)-1)/(p-1))|0);
    }
    if(reverse){
        l.reverse();
    }
    return l;
}

exports.multinomial = multinomial;
exports.binomial = binomial_gen;
exports.restricted_partitions = restricted_partitions;
exports.xi_degrees = xi_degrees;


function * WeightedIntegerVectors(n, l){
    /*
    Iterate over all ``l`` weighted integer vectors with total weight ``n``.

    INPUT:

    - ``n`` -- an integer
    - ``l`` -- the weights in weakly decreasing order

    EXAMPLES::

        sage: from sage.combinat.integer_vector_weighted import iterator_fast
        sage: list(iterator_fast(3, [2,1,1]))
        [[1, 1, 0], [1, 0, 1], [0, 3, 0], [0, 2, 1], [0, 1, 2], [0, 0, 3]]
        sage: list(iterator_fast(2, [2]))
        [[1]]

    Test that :trac:`20491` is fixed::

        sage: type(list(iterator_fast(2, [2]))[0][0])
        <type 'sage.rings.integer.Integer'&gt;
    */
    if(n < 0){
        return
    }

    if(!l){
        if(n === 0){
            yield []
        }
        return
    }
    if(l.length === 1){
        if(n % l[0] === 0){
            yield [Math.floor(n / l[0])];
        }
        return
    }

    let k = 0;
    let cur = [Math.floor(n /l[k]) + 1];
    let rem = n - cur[cur.length - 1] * l[k]; // Amount remaining
    while(cur.length > 0){
        cur[cur.length-1] --;
        rem += l[k];
        if(rem === 0){
            yield [...cur, ... Array(l.length - cur.length).fill(0)];
        } else if(cur[cur.length - 1] < 0 || rem < 0){
            rem += cur.pop() * l[k];
            k --;
        } else if(l.length === cur.length + 1){
            if(rem % l[l.length -1] === 0){
                yield [...cur, Math.floor(rem / l[l.length - 1])];
            }
        } else {
            k ++;
            cur.push(Math.floor(rem / l[k]) + 1);
            rem -= cur[cur.length - 1] * l[k];
        }
    }
}
exports.WeightedIntegerVectors = WeightedIntegerVectors;

// console.log(Array.from(WeightedIntegerVectors(3, [2,1,1])));