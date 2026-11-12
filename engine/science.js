export const complement = {A:'T', T:'A', G:'C', C:'G', N:'N'};
export const span = ([a,b]) => `[${a}, ${b})`;
export const flip = ([a,b],length) => [length-b,length-a];
export const revcomp = s => [...s].reverse().map(c=>complement[c]).join('');
export const inSpan = (i,[a,b]) => i>=a && i<b;

export function boundary(start,end,shift,L) {
  if(![start,end,shift,L].every(Number.isInteger)||start<0||start>=end||start+shift<0||L<0)
    throw new RangeError('Enter non-negative coordinates and an end greater than the start.');
  const a=start+shift,b=end+shift, width=b-a;
  const observed=Math.max(0,Math.min(b,L)-a);
  return {unshifted:[start,end], anchor:[a,b], read:[0,L], width, observed,
    outside:width-observed, contained:b<=L, minimumReadEnd:b};
}

export function compare(record, intervalA, intervalB, reverse=false) {
  const [origin,stop]=record.domain;
  for(const [a,b] of [intervalA,intervalB]) {
    if(!Number.isInteger(a)||!Number.isInteger(b)||a<origin||a>=b||b>stop)
      throw new RangeError(`Use a non-empty interval inside ${span(record.domain)}.`);
  }
  const [a,b]=intervalA, [c,d]=intervalB;
  const strA=record.segment.slice(a-origin,b-origin);
  const strB=record.segment.slice(c-origin,d-origin);
  const overlap=Math.max(0,Math.min(b,d)-Math.max(a,c));
  return {
    a:reverse?flip(intervalA,record.length):intervalA,
    b:reverse?flip(intervalB,record.length):intervalB,
    lengthA:strA.length,lengthB:strB.length,overlap,
    onlyA:strA.length-overlap,onlyB:strB.length-overlap,
    identicalString:strA===strB, identicalInterval:a===c&&b===d,
    sequenceA:reverse?revcomp(strA):strA,sequenceB:reverse?revcomp(strB):strB
  };
}

export function validate(audit,molecule) {
  const results=[];
  function assert(name,condition) {results.push({test:name,pass:Boolean(condition)});}
  const example=boundary(210,222,3,224);
  assert('Shifted edge: 11/12 bases inside',example.observed===11&&!example.contained);
  assert('One extra base restores containment',boundary(210,222,3,225).contained);
  assert('Unshifted anchor fits',boundary(210,222,0,224).contained);
  assert('Shifted inner anchor fits',boundary(203,215,3,224).contained);
  const record=audit.cases[0], second=audit.cases[1];
  const first=compare(record,record.rules[2].span,record.rules[1].span);
  assert('First read: 42 vs 87 nt, 45 exclusive',first.lengthA===42&&first.lengthB===87&&first.onlyB===45);
  const agreement=compare(second,second.rules[2].span,second.rules[1].span);
  assert('Second read: labels agree and strings differ',
    second.rules[2].native===second.rules[1].native&&!agreement.identicalString&&agreement.onlyA===33);
  const alt=compare(second,second.rules[2].span,second.rules[1].otherOptimal[0]);
  assert('Alternative optimal span restores coordinate and sequence agreement',alt.identicalInterval&&alt.identicalString);
  const reversed=compare(record,record.rules[2].span,record.rules[1].span,true);
  assert('Reverse-complement preserves overlap and round-trips',reversed.overlap===first.overlap&&
    revcomp(reversed.sequenceA)===first.sequenceA&&String(flip(reversed.a,record.length))===String(first.a));
  assert('Atomic model preserves archived 42-nt sequence',molecule.seq===first.sequenceA&&molecule.points.length===1716);
  return results;
}
