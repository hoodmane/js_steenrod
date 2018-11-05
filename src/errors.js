function invalidInput(fn, args, msg){
    return new Error(`Invalid input ${fn}(${args.map(JSON.stringify).join(", ")}): ${msg}`);
}
exports.invalidInput = invalidInput;


function inputNotNonnegativeInteger(fn, args, problem){
    return invalidInput(fn, args, `${JSON.stringify(problem)} is not a nonnegative integer`)
}
exports.inputNotNonnegativeInteger = inputNotNonnegativeInteger;

function checkNonnegativeInteger(fn, args, input_to_check){
    if(!Number.isInteger(input_to_check) || input_to_check < 0){
        throw inputNotNonnegativeInteger(fn, args, input_to_check);
    }
}
exports.checkNonnegativeInteger = checkNonnegativeInteger;
