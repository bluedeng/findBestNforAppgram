# knn_query
build for knn_query

2015-11-16: coding simple_index_knn_query, double_index_knn_query, multi_index_knn_query

2015-11-18: building the project and update to github

2015-11-19: simple_index_knn_query: ok£»
						double_index_knn_query: max_ed_average not still the same while processing the same dataset
						
2015-11-24: coding double_index_test for comparing with double_index_knn_query 
						in order to find out why exsit different max_ed_average
						
2015-11-25: important update from xiaoli
						good bye older older version, come newer version
						all works should be update because of the new one
						until now, the update not cover: simple, double, multi_index_knn_query, double_index_test
						
2015-11-28: we can use the simple_index_knn_query in the new version now!

2015-12-2:  now, we can use double_index_knn_query and multi_index_knn_query in the new version now!
						we also modify the SeqDB file for index_building initial dataCount[] 
																							and adds the old_version_knn_postprocess
															
2015-12-22: need to modify the simple_index_knn_query process and add stringQuery method.

2015-12-23: update codes between WINOS and CENTOS

2015-12-25: update codes of simple_index_knn_query for ed calculating
														SeqDB by adding function stringQuery for simple_index_knn_query
														(Do not use old_version_knn_postprocess in simple_index now.)