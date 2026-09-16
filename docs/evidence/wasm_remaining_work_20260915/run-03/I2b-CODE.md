# I2b emitted scan before timing

Artifact SHA256:
`86f85d75743d8661a5633326b2cb84862e5c986d1a26f34da41d0de9c9dfbab2`.
`work/I2b/function-indices.json` binds the actual function 1368. The group
still calls it at the original scan-eligible HSP boundary. The I2b helper has
1,313 WAT lines and 18 locals (I1: 1,264 and 20); size alone is not a benefit.

The separate normal-tier print run produced Liftoff and TurboFan code, with
byte-identical oracle output. Its time is excluded from performance evidence.
The TurboFan section begins at line 5374 of
`work/I2b-code-print/MjeNMV.MelaMJNV.tlosatx/threaded-n8/stdout.txt`.

The untraced scan computes the helper address in r8 at +0x191/+0x195.
The score-contender branch at +0x1d8 goes to +0x46f; coordinate loads at
+0x47f/+0x489/+0x49b use that retained r8. I1's second multiply/add on the
contender path is absent. The prefix truncation bounds check and saturating
subtraction remain on rejected-score jumps (+0x1de..+0x1f2), and the ordinary
per-next-helper multiply remains. This is not elimination of all bounds or
address calculations, nor a quantitative explanation of total elapsed time.

The preserved exact differential, fresh raw/thread gate and fixed normal
runtime reuse screen determine correctness and usefulness separately.
