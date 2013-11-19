var COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'};

function reverseComplement(string) {
  var result = "", complementsLocal = COMPLEMENTS;
  
  for(var i = string.length; i--;) {
    result += complementLsocal[string[i]];
  }
  
  return result;
}
