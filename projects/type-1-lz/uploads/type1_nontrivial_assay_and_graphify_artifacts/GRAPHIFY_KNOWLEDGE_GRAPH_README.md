# Graphify knowledge graph: Type-1, N=3 Landau-Zener

This folder contains a Graphify-ready knowledge graph summarizing the current state of the Type-1, N=3 Landau-Zener project.

## Files

- `graphify_type1_n3_knowledge_graph.json`: primary Graphify-ready node/edge JSON.
- `graphify_type1_n3_knowledge_graph.graphml`: GraphML version for import into graph tooling.
- `graphify_type1_n3_knowledge_graph.dot`: Graphviz DOT source with clusters.
- `graphify_type1_n3_knowledge_graph.svg`: rendered vector graph.
- `graphify_type1_n3_knowledge_graph.pdf`: rendered graph PDF.

## Node status tags

- `baseline`: core ODE/Type-1 data.
- `formalized`: algebraic/geometric object has a defined role in the project notes.
- `proved` / `proved in notes`: theorem-level or proof-record status in current internal notes.
- `validated`: numerical or bookkeeping validation.
- `target` / `proposed`: component or assay target not yet independently established.
- `falsified`: collapse or ansatz ruled out by algebraic or numerical checks.

## Critical interpretation

The graph intentionally separates:

1. the validated segmented-product bookkeeping assay;
2. the falsified point-local Airy collapse of projected Q4 windows;
3. the proposed nontrivial assay levels needed to construct projected Q4 window objects independently;
4. the selector insertion theory, which applies to active W4 selector points and should not be conflated with ordinary projected Q4 windows.
