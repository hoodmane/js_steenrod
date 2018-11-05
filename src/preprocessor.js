let esprima = require("esprima");
let estraverse = require("estraverse");
let escodegen = require("escodegen");

function getContextMemberNode(id){
    let node = {};
    node.type = "MemberExpression";
    node.object = {
        "type" : "Identifier",
        "name" : "context"
    };
    node.property = {
        "type" : "Identifier",
        "name" : id,
        "member" : true
    };
    return node;
}

let binary_operator_to_fn_map = {
    "+" : "ADD",
    "-" : "SUB",
    "*" : "MULT",
    "==": "EQUAL"
};

let unary_operator_to_fn_map = {
    "-" : "MINUS"
};

let type_to_enter_fn_map = {};

type_to_enter_fn_map.BinaryExpression = function(node) {
    // Only replace operators in operator_to_fn_map
    if (!binary_operator_to_fn_map.hasOwnProperty(node.operator)) {
        return;
    }
    node.type = "CallExpression";
    node.callee = getContextMemberNode(binary_operator_to_fn_map[node.operator]);
    node.arguments = [node.left, node.right];

    // reset unnecessary properties
    node.left = null;
    node.right = null;
    node.operator = null;

    return node;
};

type_to_enter_fn_map.UnaryExpression = function(node) {
    // Only replace operators in operator_to_fn_map
    if (!unary_operator_to_fn_map.hasOwnProperty(node.operator)) {
        return;
    }
    node.type = "CallExpression";
    node.callee = getContextMemberNode(unary_operator_to_fn_map[node.operator]);
    node.arguments = [node.argument];

    // reset unnecessary properties
    node.prefix = null;
    return node;
};

let Sq_re = /^(Sq|P|Q|bP|pst)$/;
let Sqn_re = /^(Sq|P|Q|bP)(\d+)$/;
type_to_enter_fn_map.Identifier = function(node){
    if(node.member){
        return node;
    }
    if(Sqn_re.test(node.name)){
        let op_n = Sqn_re.exec(node.name);
        let op = op_n[1];
        let n  = op_n[2];
        node.type = "CallExpression";
        node.callee = getContextMemberNode(op);
        node.arguments = [
            {
                "type": "Literal",
                "value": Number(n),
                "raw": n
            }
        ];
        node.name = null;
    } else if(Sq_re.test(node.name)){
        Object.assign(node, getContextMemberNode(node.name));
        node.identifier = null;
    } else if(node.name === "b"){
        node.type = "CallExpression";
        node.callee = getContextMemberNode("b");
        node.arguments = [];
        node.identifier = null;
        node.name = null;
    }
    return node;
};

//console.log(JSON.stringify(esprima.parse("-Sq2"),null,'\t'));
//console.log(JSON.stringify(esprima.parse("a.b.c"),null,'\t'));

function preprocess(code){
    let ast  = esprima.parse(code);
    estraverse.traverse(ast, {
        enter: function(node) {
            if(node.type === "MemberExpression"){
                node.property.member = true;
            }
            if(type_to_enter_fn_map[node.type]){
                type_to_enter_fn_map[node.type](node);
            }
        }
    });
    let modified_code = escodegen.generate(ast);
    return modified_code;
}

//preprocess("range(5).map(Sq)");
//preprocess("range(5).map(context.Sq)");
//console.log(JSON.stringify(esprima.parse("b(c)"),null,'\t'));
//console.log(preprocess("-Sq(2)"));

module.exports = preprocess;
