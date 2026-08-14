# Reaction connectivity: NFsim changes for mRNA-like molecules

## Purpose

- This note explains the changes that were made to NFsim so that it can simulate
  mRNA-like polymer molecules, such as an mRNA carrying many translating ribosomes.
- The work was started in the [rasilab/nfsim](https://github.com/rasilab/nfsim) fork and
  later merged into [RuleWorld/nfsim](https://github.com/RuleWorld/nfsim).
- The note describes the code as it stands on the current `master`, and flags where later
  upstream commits changed the behavior of the merged feature.

## How the work reached the upstream repository

- Pull request #30 (branch `connectivity`), merge commit `a7d5377`, merged 2021-07-01 by
  Ali Sinan Saglam. 29 files, 2079 insertions and 223 deletions. This is the substantive
  change.
- Pull request #32 (branch `rasilab-cleanup_output`), merge commit `ccef9bf`, merged
  2022-01-13. Only [src/NFsim.cpp](../src/NFsim.cpp) and
  [src/NFinput/NFinput.cpp](../src/NFinput/NFinput.cpp), 59 insertions and 50 deletions.
  Command line tidy-up.
- Original development ran from October 2017 to December 2018 in the fork. The fork
  branches carry a `README.org` with dated design notes by the author; that file was
  removed before the merge.

## The problem being solved

- An mRNA with many bound ribosomes is a single large connected complex.
- After any rule fires, stock NFsim asks every molecule in the product complex to re-test
  its membership in every reaction its `MoleculeType` participates in.
- In a TASEP-style translation model the ribosome and the mRNA each appear in more than a
  thousand rules, while typically only one to five propensities actually change.
- The per-firing cost is therefore approximately
  `|products| x |rules touching that MoleculeType| x cost(pattern match)`, where the
  pattern match is a recursive bond traversal rather than a constant-time operation.
- The 2018-08-01 design note states the diagnosis plainly: the neighborhood search was not
  the bottleneck, the membership update was.

## Two approaches were tried, and only one shipped

- **Polymer neighborhood search**, developed July to August 2018. A `MoleculeType` could
  be tagged as polymeric, and the bonded-neighbor breadth search stopped at a bounded
  distance so that the product list did not grow with mRNA length. The author measured
  approximately a factor of two and judged it not worth the complexity. A series of
  commits titled "Remove polymer ... crud" (`b3aaf78`, `391c8c7`, `6596428`, `f66e727`,
  among others, December 2018) stripped it out before the merge. The current `master`
  contains no polymer code.
- **Static reaction connectivity inference**, developed from August 2018 onward. This is
  what shipped. The author measured close to a hundredfold speedup on the translation
  models.

## How the shipped mechanism works

### Parse phase

- `NFinput::initReactionRules` in [src/NFinput/NFinput.cpp](../src/NFinput/NFinput.cpp)
  now retains every reactant and product `TemplateMolecule` for each rule, rather than
  only the first per connected component.
- A new parser, `NFinput::readTemplatePattern` (defined at line 4158), builds
  `TemplateMolecule` objects for **product** patterns. Products are normally realized by
  transformations and never turned into templates, so this parser exists purely to give
  the connectivity analyzer a structural description of what each rule produces.
- The BNGL XML `<Map>` block, holding `<MapItem sourceID=... targetID=.../>` entries, is
  read at lines 1586 to 1613. Each reactant template is paired with its product
  counterpart through `setMappedPartner`, so that every reactant template knows the
  template describing its own post-firing state.
- The two lists are handed to the rule by
  `ReactionClass::setAllReactantAndProductTemplates` in
  [src/NFcore/reactionClass.cpp](../src/NFcore/reactionClass.cpp) at line 294, and stored
  as `allReactantTemplates` and `allProductTemplates` in
  [src/NFcore/NFcore.hh](../src/NFcore/NFcore.hh) at lines 1415 and 1416.
- Every `Transformation` in
  [src/NFreactions/transformations/transformation.hh](../src/NFreactions/transformations/transformation.hh)
  gains a `TemplateMolecule` back-pointer, exposed through the pure virtual
  `getTemplateMolecule()` at line 218. Semantically, the transformation is the operation
  the rule performs and the template is the pattern context in which it is performed.
- A new factory method, `genUnbindingTransform2`, creates an empty transformation carrying
  the binding partner's component index and template, so that the partner side of an
  unbinding is visible to the analyzer.

### Inference phase

- Under the `-connect` switch, `System::prepareForSimulation` in
  [src/NFcore/system.cpp](../src/NFcore/system.cpp) at lines 680 to 722 allocates a dense
  boolean matrix `connectedReactions` over all rules and calls
  `ReactionClass::identifyConnectedReactions` for each rule.
- This runs before molecules are added to reactant lists. It is a purely static analysis
  of rule structure, costing order `N^2` template compatibility tests and order `N^2` bits
  of storage, performed once.
- Rule A is marked as connected to rule B when firing A could change B's reactant set. The
  relation is directional and not symmetric.
- The predicate is `ReactionClass::isReactionConnected` at
  [src/NFcore/reactionClass.cpp](../src/NFcore/reactionClass.cpp) line 254. It first calls
  `TransformationSet::checkConnection`, which walks each transformation's template, then
  that template's mapped product partner, then any added-molecule template. It then sweeps
  the rule's full reactant and product template lists.
- The structural primitive is `TemplateMolecule::isTemplateCompatible` in
  [src/NFcore/templateMolecule.cpp](../src/NFcore/templateMolecule.cpp) at line 1565. It
  answers whether a single concrete molecule could satisfy both patterns at once. Two
  templates are incompatible when they are different molecule types, sit in different
  compartments, or make contradictory claims about a shared component: occupied against
  empty, empty against bonded, conflicting component-state constraints, or explicit bonds
  to different partners. Patterns whose constrained component sets are disjoint are
  treated as compatible.
- The test deliberately over-approximates and returns true when uncertain, so the inferred
  graph is a superset of the true connectivity. The 2018-08-16 design note records this as
  a conscious tradeoff to avoid ever losing a genuine dependency.

### Simulation phase

- `ReactionClass::fire` in
  [src/NFcore/reactionClass.cpp](../src/NFcore/reactionClass.cpp) decides per product
  molecule whether the fast path is safe, at lines 571 to 579. It requires connectivity to
  be enabled, no indirect products, and the molecule to be a direct product.
- The call is routed through `Molecule::updateRxnMembership` in
  [src/NFcore/molecule.cpp](../src/NFcore/molecule.cpp) at lines 154 to 162, which
  dispatches to `MoleculeType::updateConnectedRxnMembership` in
  [src/NFcore/moleculeType.cpp](../src/NFcore/moleculeType.cpp) at lines 580 to 611.
- That loop still walks the `MoleculeType` reactions in native order, so that container
  mutation order matches the full-refresh path exactly, but it skips any rule for which
  `System::areReactionsConnected` returns false. The lookup is constant time on the
  precomputed matrix.
- The saving is the entire cost of `tryToAdd`, and therefore of the recursive
  `TemplateMolecule::compare`, for every unconnected rule.
- Note that the complex traversal itself was **not** made connectivity-aware. The breadth
  first search still returns the same molecule set; what changed is the work done per
  returned molecule.

## Command line switches

- `-connect` infers network connectivity before starting the simulation. Default is off.
  It requires no modification to BioNetGen or PySB.
- `-printconnected` writes the inferred connectivity of each reaction to a file.
- `-trackconnected` writes out which reaction rates actually changed after each firing.
- `-trackrxnnum` logs reaction numbers instead of names, to keep the log small.
- `-rxnlog [filename]` writes firing times and participating molecules as JSON. The
  expected extension is `.nfevent.json`, and the schema is in
  [nfevent.schema.json](../nfevent.schema.json).
- `-logbuffer [int]` sets how many firings to buffer between writes to the reaction log.
- `-maxcputime` sets a wall-clock ceiling in seconds. Default is no limit.
- `-utl [int]`, the universal traversal limit, predates this work but is the knob the mRNA
  test cases rely on to bound complex traversal.

## Output files

| Switch | File | Contents |
|---|---|---|
| `-rxnlog f.nfevent.json` | as given | JSON with simulation info, molecule types, initial state, and a `firings` array |
| `-rxnlog` and `-printconnected` | `<log minus 4 chars>_connectedlist.tsv` | `rxn_id rxn_name connected_rxn_id con_rxn_name`, static, requires `-connect` |
| `-rxnlog` and `-trackconnected` | `<log minus 4 chars>_connected.tsv` | `line rxn mol mol_id con_rxn old_a new_a`, per runtime propensity change |
| `-printmoltypes` | `<out>.molecule_type_list.tsv` | molecule type rows, plus `last_rxn_firing_time`, `simulated_time` and `cpu_time` |
| `-printrxncounts` | `<out>.rxn_list.tsv` | `rxn n_firings name` |

- Diffing the rules whose propensity changed against the rules marked as connected is the
  direct way to detect an under-connected or over-connected matrix.

## Validation models

- [validate/basicModels/v20.bngl](../validate/basicModels/v20.bngl) is a small enzyme
  phosphorylation model, added with pull request #30, that exercises every new switch at
  once on a system with a known ordinary differential equation answer.
- [validate/basicModels/v21.bngl](../validate/basicModels/v21.bngl) is the mRNA model. It
  declares `mrna(p0~free~blocked,p1..p9)` and `ribosome(asite)`, with rules `initiate`
  using `TotalRate`, `elongate_1` through `elongate_9`, and `terminate`. The
  `p0~free` and `p0~blocked` states encode the ribosome footprint: position `p0` frees
  only after the first elongation step, so a second ribosome cannot initiate on top of the
  first. This is a ten-codon TASEP.
- `v21.bngl` was added by Ali Sinan Saglam in `23c62c5` on 2021-07-02, one day after the
  pull request #30 merge, as part of integration rather than by the original author.
- Models `r21` through `r24` extend the mRNA model with add and delete operations, and
  `r25` and `r26` check `-rxnlog` and `-logbuffer` in isolation. All mRNA runs pass
  `-utl 3 -gml 1000000`.
- [validate/validate.py](../validate/validate.py) asserts that the NFsim deviation from
  the ordinary differential equation solution is no worse than the stochastic simulation
  algorithm deviation by more than a tolerance of 0.35.

## What upstream changed afterwards, and why it matters

- The order preservation and the broadened predicate described above are **not** the
  merged rasilab code. They come from two commits by Bill Hlavacek on 2026-05-11:
  `23436e2` "fix: preserve seeded trajectories with connected updates" and `124f8b7`
  "fix: narrow connectivity reachability for add-only paths".
- These commits trade away much of the original speedup in exchange for exact same-seed
  reproducibility.
- As merged, `ReactionClass::isReactionConnected` was a single line delegating to
  `TransformationSet::checkConnection`, and `checkConnection` pre-filtered on the specific
  changed component index through `areMoleculeTypeAndComponentPresent`. That
  component-index filter was the sharpest part of the heuristic.
- The filter is now gone. `ReactionClass::areMoleculeTypeAndComponentPresent` at
  [src/NFcore/reactionClass.cpp](../src/NFcore/reactionClass.cpp) line 675 and
  `TemplateMolecule::isMoleculeTypeAndComponentPresent` at
  [src/NFcore/templateMolecule.cpp](../src/NFcore/templateMolecule.cpp) line 2129 are
  still defined but have no callers anywhere in the tree.
- **The matrix is now denser.** `isReactionConnected` returns true whenever any reactant
  template, or for non-synthesis rules any product template, of the firing rule is
  compatible with any reactant template of the target rule, regardless of which components
  the rule actually touches. Because `isTemplateCompatible` treats disjoint constrained
  component sets as compatible, most rules mentioning a ribosome end up connected to most
  other rules mentioning a ribosome.
- **The fast path is often skipped entirely.** `ReactionClass::fire` builds a
  `directProductSet` from the mapped molecules and sets `hasIndirectProducts` if any live
  product is absent from it, at lines 485 to 503, then requires `!hasIndirectProducts`
  before taking the connected path. Because `getListOfProducts` expands breadth-first out
  to the traversal limit, any firing whose neighborhood reaches past the explicitly mapped
  molecules falls back to the full updater for every molecule. On a polymer backbone,
  which continues in both directions from any ribosome, that is the normal case rather
  than the exception.
- The guardrail for these changes is two same-seed parity tests,
  `test_connectivity_preserves_seeded_tlbr_trajectory` and
  `test_connectivity_preserves_seeded_local_function_trajectory`, in
  [validate/validate.py](../validate/validate.py) at lines 289 to 302.
- **Practical implication.** The `-connect` switch on the current `master` is correct and
  trajectory-identical, but its speedup on the polymer workload it was written for is
  likely far below the near-hundredfold measured in 2018. Anyone reviving this for mRNA
  work should first measure the `hasIndirectProducts` hit rate and the density of
  `System::connectedReactions`. Those are the two levers to relax, with the parity tests
  as the check.

## Known defects, all present on the current master

- `System::start` is declared as a public member in
  [src/NFcore/NFcore.hh](../src/NFcore/NFcore.hh) at line 522, but `System::sim()` in
  [src/NFcore/system.cpp](../src/NFcore/system.cpp) declares a local `clock_t start` at
  line 947 that shadows it and is the one actually assigned. The member is never
  initialized, yet [src/NFcore/reactionClass.cpp](../src/NFcore/reactionClass.cpp) at line
  631 computes `current_cpu_time` from `this->system->start`. The `cpu_time` value written
  into `.molecule_type_list.tsv` is therefore computed against an uninitialized value
  whenever the reaction log is enabled.
- `current_cpu_time` is refreshed only inside an `if(verbose)` block in the main loop of
  [src/NFcore/system.cpp](../src/NFcore/system.cpp) at lines 1003 to 1016, which makes
  `-maxcputime` ineffective unless `-v` or `-rxnlog` is also passed. The release notes in
  [README.md](../README.md) record the option as disabled in version 1.2.2.
- The derived filenames for `-trackconnected` and `-printconnected` in
  [src/NFsim.cpp](../src/NFsim.cpp) at lines 685 to 704 chop the last four characters of
  the log name using `std::string::replace`, which mutates the string in place. With the
  recommended `.nfevent.json` extension this strips `json`, and passing both switches
  makes the second block operate on the already mutated name, producing
  `_connected_connectedlist.tsv`.
- `TransformationSet::checkConnection` in
  [src/NFreactions/transformations/transformationSet.cpp](../src/NFreactions/transformations/transformationSet.cpp)
  returns false on the first `REMOVE` transformation it encounters, at lines 988 to 992,
  rather than continuing the loop. The fallback loops in `isReactionConnected` mask this
  today, but it is presumably the origin of the bug recorded in
  [README.md](../README.md) as "`-connect` would break with models that have remove
  operations".
- `ReactionClass::appendPreConnectedRxn` is declared in
  [src/NFcore/NFcore.hh](../src/NFcore/NFcore.hh) at line 1374 but is never defined or
  called. `ReactantList::getWriteableMappingSet` in
  [src/NFreactions/reactantLists/reactantList.cpp](../src/NFreactions/reactantLists/reactantList.cpp)
  at line 152 is defined but never called. Both are leftovers from the original branch.

## References

- Release notes for the feature and its bugfixes are in [README.md](../README.md).
- The reaction log format is specified in [nfevent.schema.json](../nfevent.schema.json).
- The shipped [NFsim manual](../NFsim_manual_v1.12.pdf) predates all of this work and does
  not describe any of these switches. They are documented only in the `-help` output and
  in the Doxygen header of [src/NFsim.cpp](../src/NFsim.cpp).
