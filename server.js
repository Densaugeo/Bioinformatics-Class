// This is just a little simple REPL for dna.js
var repl = require('repl');
var dna = require('dna.js');

console.log('Hello, welcome to dna.js\' REPL!  The library is available in the \'dna\' object, and all of its properties can be referenced directly by name from here.');

var cli = repl.start({});
cli.context.dna = dna;
for(var i in dna) cli.context[i] = dna[i];
