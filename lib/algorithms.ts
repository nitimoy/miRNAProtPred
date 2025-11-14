// Boyer-Moore Search Algorithm
export function preprocessBadCharacterTable(pattern: string): Map<string, number> {
  const table = new Map<string, number>()
  for (let i = 0; i < pattern.length; i++) {
    table.set(pattern[i], pattern.length - i - 1)
  }
  return table
}

export function boyerMooreSearch(text: string, pattern: string): number[] {
  const m = pattern.length
  const n = text.length

  if (m > n) {
    return []
  }

  const table = preprocessBadCharacterTable(pattern)
  const occurrences: number[] = []

  let i = m - 1
  while (i < n) {
    let j = m - 1
    while (j >= 0 && text[i] === pattern[j]) {
      i--
      j--
    }

    if (j === -1) {
      occurrences.push(i + 1)
    }

    const badCharSkip = table.get(text[i]) ?? m
    i += Math.max(m - j, badCharSkip)
  }

  return occurrences
}

// Genetic code dictionary
export const geneticCode: { [key: string]: string } = {
  A: 'GCU',
  R: 'CGU',
  N: 'AAU',
  D: 'GAU',
  C: 'UGU',
  Q: 'CAA',
  E: 'GAA',
  G: 'GGU',
  H: 'CAU',
  I: 'AUU',
  L: 'UUA',
  K: 'AAA',
  M: 'AUG',
  F: 'UUU',
  P: 'CCU',
  S: 'UCU',
  T: 'ACU',
  W: 'UGG',
  Y: 'UAU',
  V: 'GUU',
  '*': 'UAA',
}

export function rnaToDNA(rnaSequence: string): string {
  return rnaSequence.replace(/U/g, 'T')
}

export function proteinToRNA(proteinSequence: string): string {
  return proteinSequence
    .split('')
    .map((aa) => geneticCode[aa] || '')
    .join('')
}
