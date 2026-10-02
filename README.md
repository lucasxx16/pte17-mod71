# The size-17 symmetric PTE obstruction modulo 71

This repository records a **computational nonexistence proof** for the following finite-field problem, together with exact certificates, source code, saved enumeration outputs, and two separate AI-agent reviews.

**Result.** If a multiset $A$ of total size $17$ in $\mathbb F_{71}$ satisfies

$$\sum_{a\in A}a^j=0\qquad(j=1,3,5,7,9,11,13,15),$$

then $A=-A$ as multisets. **Repeated elements are allowed and counted.** Thus no nontrivial multiset of the specified kind exists.

中文：本仓库给出“模71的17元多重集八个奇矩同时为零必取负对称”的计算辅助证明。完整中文证明见 [D 的证明](proof/round1/D/report.txt)，独立审稿见 [C](proof/round1/C/review_D_global.txt) 和 [E](proof/round1/E/review_D_global.txt)。下述复跑入口会在独立副本中执行，不改动保存的证据。

## Proof structure

Newton's identities give a monic root polynomial $P(X)=f(X)-c$, where $f$ is odd of degree 17. A nontrivial solution would have $c\ne0$, no zero root, and no pair of opposite roots. If it has $r$ distinct roots, evenness of $f'$ implies $2(17-r)\le16$, so $r\ge9$.

The proof partitions all possible multiplicities into two exhaustive branches:

| Branch | Complete search | Saved outcome |
| --- | --- | --- |
| A root has multiplicity at least 4 | Scale that root to 1 and enumerate $P(X)=17\int_1^X(T^2-1)^3R(T^2)\,dT$, with $R$ monic of degree 5. All $71^5=1,804,229,351$ parameter choices are covered through $71^4$ histograms. | 52,910 screen survivors; 6,534 with nonzero constant; none splits completely, counting multiplicity. |
| All root multiplicities are at most 3 | The signed multiplicity vector $v\in\mathbb Z^{35}$ lies in the full lattice $L=\{v:Hv=0\pmod{71}\}$, with $\|v\|_1=17$ and $|v_i|\le3$. Since $m^2\le4m-3$ for $m=1,2,3$ and $r\ge9$, $\|v\|_2^2\le41$. Enumerate this covering region using certified integer bounds. | 15,207,867 leaves; 7,641,386 odd-parity leaves; zero target vectors. |

Here $H_{j,a}=a^j$ for $j\in\{1,3,\ldots,15\}$ and $a=1,\ldots,35$. An exact certificate proves that the saved basis generates **all** of $L$, of determinant $71^8$. Rational Gram–Schmidt checks certify the fixed-point constants and conservative integer pruning. Floating-point LLL was only used to discover a convenient basis; it is not a correctness dependency of reproduction.

The two branches exclude every nontrivial multiset. The resulting local obstruction implies that the constant polynomial difference of an integer size-17 symmetric ideal Prouhet–Tarry–Escott pair $(A,-A)$ is divisible by 71. This is not a classification of integer solutions.

The problem context is [arXiv:2304.11254, §2.2](https://arxiv.org/html/2304.11254v1#S2.SS2). No priority claim is made.

## Reproduce

The recorded environment was **Windows, Python 3.14.7, NumPy 2.5.2, and PowerShell 7 with C#/.NET `Add-Type`**. The full audit requires PowerShell 7+ (`pwsh`); Windows PowerShell 5.1 lacks the audit's `Add-Type -CompilerOptions` parameter. No network access is required after dependencies are installed. A reproduction uses one process at a time and usually takes several minutes; runtime depends on the machine.

From the repository root:

```powershell
python -m pip install -r requirements.txt
python -E scripts/reproduce.py
```

The default run checks file integrity, copies the evidence into a new ignored `.runs/` directory, verifies all three exact lattice certifications, repeats the full high-multiplicity and lattice searches, runs the different high-multiplicity implementation and other independent checks, and compares the full semantic outputs with the saved evidence. Logs and a machine-readable `run_summary.json` remain in that run directory. Expected final status: `PASS`, with `full_enumerations: true`.

The public package passed this complete 12-step reproduction on 2026-10-02: [machine-readable validation record](validation/public_reproduction.json). Both full search outputs matched the frozen evidence, including all lattice layer counts.

For a short certificate/control check without the exhaustive enumerations:

```powershell
python -E scripts/reproduce.py --quick
```

A quick pass **does not rerun the nonexistence searches**. It checks the saved exact lattice certificates, a known short-vector positive control, the integer-square-root audit, the reference witness checker's controls, and file integrity.

To select a PowerShell executable explicitly:

```powershell
python -E scripts/reproduce.py --powershell pwsh.exe
```

Python assertions are part of the auditors. Do not use `-O`; the wrapper invokes child Python processes with `-E` and without optimization. To check only the saved proof-file hashes, run `python -E scripts/check_integrity.py`.

For manual reproduction, copy the entire `proof/` directory first. In the copy, run `verify_lattice_certificate.py`, `check_multiplicity4.py`, and `run_integer_lattice.ps1` from `round1/D/`. Keep the relative `round1/D`, `round1/C/audit_D`, and `round1/E/audit_D` layout intact. Original scripts write outputs in their working directories.

## Evidence and scope of review

| Evidence | Location |
| --- | --- |
| Complete proof, primary source, basis certificate, and original outputs | [proof/round1/D/](proof/round1/D/) |
| Adversarial review; independent Gram-matrix, Hasse-multiplicity and positive-path checks; saved full replays | [proof/round1/C/review_D_global.txt](proof/round1/C/review_D_global.txt), [audit code](proof/round1/C/audit_D/) |
| Separate review; different full high-multiplicity enumeration; Bareiss determinant and explicit rational Gram–Schmidt checks; overflow bounds; exhaustive square-root checks; saved full replays | [proof/round1/E/review_D_global.txt](proof/round1/E/review_D_global.txt), [audit code](proof/round1/E/audit_D/) |
| Original multiset witness verifier and small control tests | [proof/reference/verify_local.py](proof/reference/verify_local.py) |
| Supplementary structural reductions, including the odd-multiplicity-class criterion and orbit counts | [E's structural report](proof/round1/E/report.txt), [E7 review](proof/round1/C/review_E7.txt), [E8 review](proof/round1/C/review_E8.txt) |

The proof was developed by an independent **GPT-6 Astra xhigh** agent (D), and separately reviewed by **GPT-6 Sol** (C) and another **GPT-6 Astra xhigh** agent (E). D's first draft was frozen before the cross-reviews; C and E did not share their initial review conclusions. The coordinator checked the decisive reports, code and outputs. These are **AI-agent reviews, not external human peer review**.

The high-multiplicity branch has a different complete implementation by E. Both reviewers also repeated the full lattice enumeration, but those lattice replays use the **same audited core enumerator**. There is no second independently programmed full lattice enumerator and no proof-assistant formalization. Correctness depends on the stated mathematics, the programs, and their execution environment.

Original proof files are retained byte-for-byte, including historical pending-review wording and elapsed times. Later C/E reviews resolve the primary draft's pending review. Supplementary reports contain historical or pending claims; they are not dependencies of the two-branch final proof. Some historical references point to exploratory files outside this public proof package. The included `lattice_records.json` supplies only the mass-18 positive control; it is not a mass-17 witness.

`proof/SHA256.json` records the exact included evidence files. Hashes check transport integrity, not mathematical correctness. The `.gitattributes` file disables newline conversion for this evidence so the hashes remain portable.
