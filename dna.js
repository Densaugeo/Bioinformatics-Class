var COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'};

// Reverse and complement DNA string
function reverseComplement(string) {
  var result = '', complements = COMPLEMENTS;
  
  // Step backwards through array for reverse iteration
  for(var i = string.length; i--;) {
    result += complements[string[i]];
  }
  
  return result;
}

// Find all matches of pattern in sequence
function findMatches(sequence, pattern) {
  var result = [];
  
  for(var last = sequence.indexOf(pattern); last !== -1; last = sequence.indexOf(pattern, last + 1)) {
    result.push(last);
  }
  
  return result;
}

// Find near matches of pattern in sequence, with no more than d mismatched nucleotides
function findApproximateMatches(sequence, pattern, d) {
  var result = [], subsequence = '', k = pattern.length, differences = 0;
  
  for(var position = 0, end = sequence.length - pattern.length; position <= end; ++position) {
    subsequence = sequence.substr(position, k);
    
    differences = 0;
    for(var j = 0; j < k; ++j) {
      if(subsequence[j] !== pattern[j]) ++differences;
    }
    
    if(differences <= d) result.push(position);
  }
  
  return result;
}

// Find all clumps of kmers which appear at least t times in any substring of length L in the given genome
function findClumps (genome, k, L, t) {
  var results = [], clumps = {}, kmers = {}, newString = '' oldString = '';
  
  for(var position = 0, end = genome.length - k; position <= end; ++position) {
    // Catalogue each kmer in kmers object, and remove old kmers from before the current substring of length L
    newString = genome.substr(position, k);
    kmers[newString] = (kmers[newString] || 0) + 1;
    if(position >= L) --kmers[genome.substr(position - L, k)];
    
    // Check for clumps
    if(kmers[newString] >= t) clumps[newString] = {};
  }
  
  for(var i in clumps) results.push(i);
  return results;
}

// Find indices of minimum G - C skew (0 based indexing)
function minimumSkew(sequence) {
  var lowestPoints = [-1], lowestValue = 0, skew = 0, skewValue = {'A': 0, 'T': 0, 'C': -1, 'G': 1};
  
  for(var i = 0, length = sequence.length; i < length; ++i) {
    skew += skewValue[sequence[i]];
    
    if(skew === lowestValue) {
      lowestPoints.push(i);
    }
    if(skew < lowestValue) {
      lowestPoints = [i];
      lowestValue = skew;
    }
  }

  return {points: lowestPoints, value: lowestValue};
}

// Find most common kmer(s) in sequence
function frequentWords(sequence, k) {
  var kmers = {}, frequentKmers = [], highestFrequency = 0, subsequence = '';
  
  for(var i = 0, end = sequence.length - k; i <= end; ++i) {
    subsequence = sequence.substr(i, k);
    
    kmers[subsequence] = (kmers[subsequence] || 0) + 1;
    if(kmers[subsequence] === highestFrequency) {
      frequentKmers.push(subsequence);
    }
    if(kmers[subsequence] > highestFrequency) {
      frequentKmers = [subsequence];
      highestFrequency = kmers[subsequence];
    }
  }
  
  return {frequency: highestFrequency, kmers: frequentKmers};
}

// Find most common kmer(s) in sequence, coutning occurences of a kmer with no more than d differences (rather than exact occurences)
/*function frequentApproximateWords(sequence, k, d) {
  var kmers = {}, frequentKmers = [], highestFrequency = 0, subsequence = '';
  
  for(var i = 0, end = sequence.length - k; i <= end; ++i) {
    subsequence = sequence.substr(i, k);
    
    kmers[subsequence] = (kmers[subsequence] || 0) + 1;
    if(kmers[subsequence] === highestFrequency) {
      frequentKmers.push(subsequence);
    }
    if(kmers[subsequence] > highestFrequency) {
      frequentKmers = [subsequence];
      highestFrequency = kmers[subsequence];
    }
  }
  
  return {frequency: highestFrequency, kmers: frequentKmers};
}*/
