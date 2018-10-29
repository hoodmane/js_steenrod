let SerreCartanBasis = require("./SerreCartanBasis");
let MilnorBasis = require("./MilnorBasis");
let Vector = require("./vector");

const jsep = require("jsep");
//jsep.addBinaryOp("^", 10);
jsep.addBinaryOp("=",1);

let binop_map = {};
let unop_map = {};
let function_map = {};
let identifier_map = {};

binop_map["*"] = function(u, v, context){
    if(u.constructor === Number){
        if(v.constructor === Number){
            return u*v;
        } else {
            return v.scale(u);
        }
    }
    return u.mult(v);
};

binop_map["+"] = function(u, v, context){
    return u.add(v);
};

function_map["Sq"] = function(args, context){
    return context.basis.Sq(args, context.p);
};

function_map["P"] = function(args, context){
    return context.basis.P(args, context.p);
};

function_map["Q"] = function(args, context){
    return context.basis.Q(args, context.p);
};

function_map["bP"] = function(args, context){
    return context.basis.bP(args[0], context.p);
};

identifier_map["b"] = function(context){
    return context.basis.b(context.p);
};

function handle_assignment(lhs, rhs, context) {
    if(!context.identifiers){
        context.identifiers = {};
    }
    if(!context.identifiers[context.p]){
        context.identifiers[context.p] = {};
    }
    if(lhs.type === "Identifier"){
        let result = evaluate(rhs, context);
        context.identifiers[context.p][lhs.name] = result;
        return result;
    } else {
        throw new Error(`Invalid left hand side of assignment. ${tree_to_string(lhs)} is not a valid identifier.`);
    }
}

function evaluate_binop(operator, left, right, context) {
    if(operator === "="){
        handle_assignment(left, evaluate(right,context), context);
    }
    return binop_map[operator](evaluate(left, context), evaluate(right, context), context);
}

function evaluate_unop(operator, argument, prefix, context) {
    return unop_map[operator](argument, prefix, context);
}

function evaluate_logical(operator, left, right, context) {
    return binop_map[operator](left, right, context);
}

function evaluate_call(callee, arguments, context) {
    return function_map[callee.name](arguments, context);
}

function evaluate_member(object, property, context) {
    throw new Error("Not implemented");
    return undefined;
}

function evaluate_array(elements, context) {
    throw new Error("Not implemented");
    return undefined;
}

function evaluate_identifier(name, context) {
    if(identifier_map.hasOwnProperty(name)){
        return identifier_map[name](context);
    } else if(context.identifiers && context.identifiers.hasOwnProperty(name)){
        return context.identifiers[name];
    }
    throw new Error("Unknown identifier");
}

function evaluate_literal(value, raw) {
    return value;
}

function evaluate(tree, context){
    switch(tree.type){
        case 'BinaryExpression' :
            // Don't recurse with evaluate in case the operator is "=" (I guess I should change the other cases like this too...)
            return evaluate_binop(tree.operator, tree.left, tree.right, context);
        case 'UnaryExpression':
            return evaluate_unop(tree.operator, evaluate(tree.argument, context), tree.prefix, context);
        case 'LogicalExpression' :
            return evaluate_logical(tree.operator, evaluate(tree.left, context), evaluate(tree.right, context), context);
        case 'CallExpression':
            return evaluate_call(tree.callee, tree.arguments.map((t) => evaluate(t, context)), context);
        case 'MemberExpression':
            return evaluate_member(evaluate(tree.object, context), evaluate(tree.property, context), context);
        case 'ArrayExpression':
            return evaluate_array(tree.elements.map((t) => evaluate(t, context)), context);
        case 'Identifier':
            return evaluate_identifier(tree.name, context);
        case 'Literal':
            return evaluate_literal(tree.value, tree.raw, context)
    }
}

function tree_to_string(tree, context){
    switch(tree.type){
        case 'BinaryExpression' :
            // Don't recurse with evaluate in case the operator is "=" (I guess I should change the other cases like this too...)
            return tree_to_string(tree.left) + " " + tree.operator + " " + tree_to_string(tree.right);
        case 'UnaryExpression':
            let arg = tree_to_string(tree.argument);
            return tree.prefix ? tree.operator + arg : arg + tree.operator;
        case 'LogicalExpression' :
            return tree_to_string(tree.left) + " " + tree.operator + " " + tree_to_string(tree.right);
        case 'CallExpression':
            return tree.callee.name + "(" + tree.arguments.map(tree_to_string).join(", ") + ")";
        case 'MemberExpression':
            return tree_to_string(tree.object) + "[" + tree_to_string(tree.property) + "]";
        case 'ArrayExpression':
            return "[" + tree.elements.map(tree_to_string).join(", ") + "]";
        case 'Identifier':
            return tree.name;
        case 'Literal':
            return tree.raw;
    }
}

function steenrod_simplify(expr, context){
    return evaluate(jsep(expr),context);
}

context = {p:2, basis : SerreCartanBasis};
console.log(steenrod_simplify("Sq(2)*Sq(2)",context));
//console.log(steenrod_simplify("x*P(1)",context));
//console.log(context);

exports.steenrod_simplify = steenrod_simplify;