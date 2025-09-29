xlog2x(x) = xlogx(x) / log(oftype(x,2))

function successes_tuple_str(yes::Int, no::Int)
    p_succ = yes / (yes + no)
    err = sqrt(p_succ * (1 - p_succ) / (yes + no))
    perc = round(100p_succ; digits=1)
    perc_err = round(100err; digits=1)
    return "conclusive: $(yes + no) | responsive: $yes ($perc ± $perc_err %) | not responsive: $no ($(round(100-perc; digits=1)) ± $perc_err %)"
end
