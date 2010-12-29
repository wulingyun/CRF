decode.exact <- function(crf)
	.Call("Decode_Exact", crf)

decode.chain <- function(crf)
	.Call("Decode_Chain", crf)

decode.tree <- function(crf)
	.Call("Decode_Tree", crf)
