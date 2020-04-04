# ComplexNetworkMedicine       
## Project updated 4th April 2020

Diseases are often caused by defective proteins, these proteins rarely operate in isolation and may have several roles in the cell. Thus over time a defective protein may be involved in several disorders, either directly or indirectly. The multiple roles leads to the concept of a disease module or cluster. This work describes how we generate overlapping clusters from complex networks to explore the dynamic nature of diseases, the genes implicated with them and the drugs used to treat them. Link clustering is vital for community detection as it enables the integration of disparate sources of data and provides a better understanding of community  hierarchy and community dynamics than non-link methods.  Furthermore,  we view not just the genes directly shared between diseases but also indirectly connected genes in the network neighborhood.  We use data and information from the STITCH protein and drug interaction databases, OMIM disease database, lists of diseases categorized by MeSH and the drugbank repository. The Gene Ontology, Disease Ontology and KEGG  provide biological validity for the disease communities. We demonstrate how the detection of overlapping clusters enables the identification of biologically plausible communities consisting of cooperating proteins. We verify their role in disease with respect to targeting drugs more effectively with expert opinion.  We have been able to identify various modules that make sense from a biological and medical perspective and validate drug repositioning candidates with clinicaltrials.gov.

![Alt Text](https://user-images.githubusercontent.com/11558110/29874468-c4ac615a-8d8e-11e7-8098-b3f18460bdf6.jpg)

If you find the data and R source code useful, please cite our paper:

K. McGarry, D. Nelson and M. Ashton, A method to explore the connectivity patterns of proteins and drugs for identifying disease communities, Accepted for publication in Springer Nature Computer Science Journal, 2020.

## References
> Barabasi, A., Gulbahce, N. and Loscalzo, J. [2011], Network medicine: a network-based approach to
human disease., Nat Rev Genet 12, 56-68.

> Barrenas, F., Chavali, S., Holme, P., Mobini, R. and Benson, M. [2009], Network properties of complex
human disease genes identified through genome-wide association studies, PLoS ONE 4(11), e8090.

> Goh, K.-I., Cusick, M. E., Valle, D., Childs, B., Vidal, M. and Barabsi, A.-L. [2007], The human disease
network, Proceedings of the National Academy of Sciences 104(21), 8685-8690.
URL: http://www.pnas.org/content/104/21/8685.abstract

> McGarry, K and Daniel, U., [2014] Computational Techniques for Identifying Networks of Interrelated Diseases, The 14th UK Workshop on Computational Intelligence, UKCI-2014, Bradford, Uk, 8th-10th Sept.

> Menche, J., Sharma, A., Kitsak, M., Ghiassian, S., Vidal, M., Loscalzo, J. and Barabasi, A. [2015], Un-
covering disease-disease relationships through the incomplete human interactome, Science 347(6224).

> Vidal, M., Cusick, M. and Barabasi, A. [2011], Interactome networks and human disease, Cell
144(6), 986-998.
