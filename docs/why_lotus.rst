What's New at Lotus?
====================

Lotus is a single-cell RNA-seq analysis pipeline, and it aims to achieve three major goals.

The first goal is to design statistically rigorous methods and robust algorithms. In this direction, we introduced a new unsupervised learning framework based on the layered structure of data points. The high-level insight is that not all cells are in the same cell state. For example, some cells may be in a stable state, while others may be intermediate cells or even low-quality cells. The purpose of the layered structure is to capture these representative cells, which we refer to as cores. Through our studies, we observe that a wide range of unsupervised algorithms, including clustering methods and visualization techniques, become more separable and more interpretable when applied to cores.

The second goal of Lotus is to make the framework fully accessible, even to users who do not know any Python or R. To achieve this, we are developing a web-based interactive projector, allowing anyone to use the tool simply by clicking buttons. For this component, we greatly welcome comments and suggestions, especially from biologists. As the next step, we plan to build an AI agent to help users operate the package more effectively.

The third goal is to integrate Lotus into multi-omics analyses. For example, our team is also working on RNA splicing analysis (https://github.com/CrossOmics/Krinon). In this direction, we see many exciting opportunities, such as single-cell level splicing analysis and other cross-omics applications. This is a long-term goal, as it requires a deep understanding of both biology and computer science. We are very excited to learn new things in both areas.
