let evaluator = require("./evaluator");
window.SerreCartanBasis = require("./SerreCartanBasis");
window.MilnorBasis = require("./MilnorBasis");
window.repl = require("repl");
window.steenrod_simplify = evaluator.steenrod_simplify;



var replServer = repl.start({
    prompt: "my-app > ",
});

