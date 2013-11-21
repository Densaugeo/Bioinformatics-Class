// Reverse and complement DNA string
function reverseComplement(string) {
  var result = '', complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'};
  
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

var changeBase = {'A': {'T': 1, 'C': 2, 'G': 3},
                  'T': {'C': 1, 'G': 2, 'A': 3},
                  'C': {'G': 1, 'A': 2, 'T': 3},
                  'G': {'A': 1, 'T': 2, 'C': 3}};

// Find most common kmer(s) in sequence, counting occurences of a kmer with no more than d differences (rather than exact occurences)
// If complements option is true, count kmers together with their complements
function frequentApproximateWords(sequence, k, d, complements) {
  if(d > 3) throw new Error('Sorry, I can only do d values of up to 3');
  console.warn('This function is complicated and untested - use at your own (high) risk')
  
  var kmers = {}, frequentKmers = [], highestFrequency = 0, subsequence = '', nearSubsequece = '', frequency = 0;
  
  for(var i = 0, end = sequence.length - k; i <= end; ++i) {
    subsequence = sequence.substr(i, k);
    
    // Make note of the subsequence, along with all similar subsequences
    kmers[subsequence] = (kmers[subsequence] || 0) + 1;
    if(d >= 1) for(var i2 = 0; i2 < k; ++i2) for(var i3 in changeBase[subsequence[i2]]) {
      nearSubsequence = subsequence.substring(0, i2) + i3 + subsequence.substring(i2 + 1);
      kmers[nearSubsequence] = (kmers[nearSubsequence] || 0) + 1;
    }
    if(d >= 2) for(var i2 = 0; i2 < k; ++i2) for(var i3 = i2 + 1; i3 < k; ++i3) {
      for(var i4 in changeBase[subsequence[i2]]) for(var i5 in changeBase[subsequence[i3]]) {
        nearSubsequence = subsequence.substring(0, i2) + i4 + subsequence.substring(i2 + 1, i3) + i5 + subsequence.substring(i3 + 1);
        kmers[nearSubsequence] = (kmers[nearSubsequence] || 0) + 1;
      }
    }
    if(d >= 3) for(var i2 = 0; i2 < k; ++i2) for(var i3 = i2 + 1; i3 < k; ++i3) for(var i4 = i3 + 1; i4 < k; ++i4) {
      for(var i5 in changeBase[subsequence[i2]]) for(var i6 in changeBase[subsequence[i3]]) for(var i7 in changeBase[subsequence[i4]]) {
        nearSubsequence = subsequence.substring(0, i2) + i5 + subsequence.substring(i2 + 1, i3) + i6 + subsequence.substring(i3 + 1, i4) + i7 + subsequence.substring(i4 + 1);
        kmers[nearSubsequence] = (kmers[nearSubsequence] || 0) + 1;
      }
    }
  }
  
  for(var i in kmers) {
    frequency = kmers[i] + (complements ? kmers[reverseComplement(i)] : 0);
    
    if(frequency === highestFrequency) {
      frequentKmers.push(i);
    }
    if(frequency > highestFrequency) {
      frequentKmers = [i];
      highestFrequency = frequency;
    }
  }
  
  return {frequency: highestFrequency, kmers: frequentKmers};
}

// Lookup table for converting DNA nucleotides to RNA
var DNA_TO_RNA = {'A': 'A', 'T': 'U', 'C': 'C', 'G': 'G'};

// Lookup table for converting RNA codons to one-letter amino acid codes
var RNA_TO_AMINES = {
  AAA: 'K', AAU: 'N', AAC: 'N', AAG: 'K',
  AUA: 'I', AUU: 'I', AUC: 'I', AUG: 'M',
  ACA: 'T', ACU: 'T', ACC: 'T', ACG: 'T',
  AGA: 'R', AGU: 'S', AGC: 'S', AGG: 'R',
  UAA: ' ', UAU: 'Y', UAC: 'Y', UAG: ' ',
  UUA: 'L', UUU: 'F', UUC: 'F', UUG: 'L',
  UCA: 'S', UCU: 'S', UCC: 'S', UCG: 'S',
  UGA: ' ', UGU: 'C', UGC: 'C', UGG: 'W',
  CAA: 'Q', CAU: 'H', CAC: 'H', CAG: 'Q',
  CUA: 'L', CUU: 'L', CUC: 'L', CUG: 'L',
  CCA: 'P', CCU: 'P', CCC: 'P', CCG: 'P',
  CGA: 'R', CGU: 'R', CGC: 'R', CGG: 'R',
  GAA: 'E', GAU: 'D', GAC: 'D', GAG: 'E',
  GUA: 'V', GUU: 'V', GUC: 'V', GUG: 'V',
  GCA: 'A', GCU: 'A', GCC: 'A', GCG: 'A',
  GGA: 'G', GGU: 'G', GGC: 'G', GGG: 'G'
}

// Lookup table for converting DNA codons to one-letter amino acid codes
var DNA_TO_AMINES = {
  AAA: 'K', AAT: 'N', AAC: 'N', AAG: 'K',
  ATA: 'I', ATT: 'I', ATC: 'I', ATG: 'M',
  ACA: 'T', ACT: 'T', ACC: 'T', ACG: 'T',
  AGA: 'R', AGT: 'S', AGC: 'S', AGG: 'R',
  TAA: ' ', TAT: 'Y', TAC: 'Y', TAG: ' ',
  TTA: 'L', TTT: 'F', TTC: 'F', TTG: 'L',
  TCA: 'S', TCT: 'S', TCC: 'S', TCG: 'S',
  TGA: ' ', TGT: 'C', TGC: 'C', TGG: 'W',
  CAA: 'Q', CAT: 'H', CAC: 'H', CAG: 'Q',
  CTA: 'L', CTT: 'L', CTC: 'L', CTG: 'L',
  CCA: 'P', CCT: 'P', CCC: 'P', CCG: 'P',
  CGA: 'R', CGT: 'R', CGC: 'R', CGG: 'R',
  GAA: 'E', GAT: 'D', GAC: 'D', GAG: 'E',
  GTA: 'V', GTT: 'V', GTC: 'V', GTG: 'V',
  GCA: 'A', GCT: 'A', GCC: 'A', GCG: 'A',
  GGA: 'G', GGT: 'G', GGC: 'G', GGG: 'G'
}

// Lookup table for converting one-letter amino acid codes to DNA codons
var AMINES_TO_DNA = {
  K: ['AAA', 'AAG'],
  N: ['AAC', 'AAU'],
  T: ['ACA', 'ACC', 'ACG', 'ACU'],
  R: ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'],
  S: ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'],
  I: ['AUA', 'AUC', 'AUU'],
  M: ['AUG'],
  Q: ['CAA', 'CAG'],
  H: ['CAC', 'CAU'],
  P: ['CCA', 'CCC', 'CCG', 'CCU'],
  L: ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
  E: ['GAA', 'GAG'],
  D: ['GAC', 'GAU'],
  A: ['GCA', 'GCC', 'GCG', 'GCU'],
  G: ['GGA', 'GGC', 'GGG', 'GGU'],
  V: ['GUA', 'GUC', 'GUG', 'GUU'],
  Y: ['UAC', 'UAU'],
  C: ['UGC', 'UGU'],
  W: ['UGG'],
  F: ['UUC', 'UUU'],
  ' ': ['UAA', 'UAG', 'UGA']
}

// Converts a given DNA string to an amine string.  Assumes sequence begins with the first nucleotide given
function dnaToAmines(sequence) {
  var result = '';
  
  for(var i = 0, end = Math.floor(sequence.length/3); i < end; ++i) {
    result += DNA_TO_AMINES[sequence.substr(3*i, 3)];
    
    if(result[i] === ' ') return result.substring(0, i);
  }
  
  return result;
}

// Searches a DNA sequence for subsequences encoding the given peptide.  Also checks reverse complement
function findPeptideInDna(sequence, peptide) {
  var results = [];
  var dnaLength = 3*peptide.length;
  var sequenceRC = reverseComplement(sequence);
  
  for(var i = 0, end = sequence.length - dnaLength; i <= end; ++i) {
    if(peptide === dnaToAmines(sequence.substr(i, dnaLength))) {
      results.push(sequence.substr(i, dnaLength));
    }
    if(peptide === dnaToAmines(sequenceRC.substr(i, dnaLength))) {
      results.push(reverseComplement(sequenceRC.substr(i, dnaLength)));
    }
  }
  
  return results;
}

// Lookup table for converting one- and three-letter amino acid codes to monoisotopic masses
var AMINES_TO_MASSES = {
  A:  71, Ala:  71,
  R: 156, Arg: 156,
  N: 114, Asn: 114,
  D: 115, Asp: 115,
  C: 103, Cys: 103,
  E: 129, Glu: 129,
  Q: 128, Gln: 128,
  G:  57, Gly:  57,
  H: 137, His: 137,
  I: 113, Ile: 113,
  L: 113, Leu: 113,
  K: 128, Lys: 128,
  M: 131, Met: 131,
  F: 147, Phe: 147,
  P:  97, Pro:  97,
  S:  87, Ser:  87,
  T: 101, Thr: 101,
  W: 186, Trp: 186,
  Y: 163, Tyr: 163,
  V:  99, Val:  99
}

// Peptide string -> monoisotopic mass
function peptideToMass(peptide) {
  var result = 0;
  
  for(var i = 0, end = peptide.length; i < end; ++i) {
    result += AMINES_TO_MASSES[peptide[i]];
  }
  
  return result;
}

// Generate full set of possible subpeptides (including the empty and original peptides)
function peptideToSubpeptides(peptide, cyclic) {
  var results = ['', peptide];
  
  // Need double peptide string for handling cyclic peptides, and save the length
  var length = peptide.length;
  peptide += peptide;
  
  for(var i = 0; i < length; ++i) for(var j = 1, end = cyclic ? length : length - i; j < end; ++j) {
    results.push(peptide.substr(i, j));
  }
  
  return results;
}

// Generate a theoretical mass spec spectrum from the given peptide
function peptideToSpectrum(peptide, cyclic) {
  var results = peptideToSubpeptides(peptide, cyclic);
  
  for(var i = 0, length = results.length; i < length; ++i) {
    results[i] = peptideToMass(results[i]);
  }
  
  results.sort(function(a, b) {
    if (a < b)
      return -1;
    if (a > b)
      return 1;
    return 0;
  });
  
  return results;
}
