

rule delta_tajima_d:
    input:
        vcf=config["vcf_file"]
    output:
        scores="results/delta_tajima_d/delta_tajima_d.scores.txt",
    params:
        window_size=config["window_size"],
        step_size=config["step_size"]
    log:
        "logs/delta_tajima_d/delta_tajima_d.log"
    script:
        "../scripts/delta_tajima_d.py"


rule tajima_d:
    input:
        vcf=config["vcf_file"]
    output:
        scores="results/tajima_d/tajima_d.scores.txt"
    params:
        window_size=config["window_size"],
        step_size=config["step_size"]
    log:
        "logs/tajima_d/tajima_d.log"
    script:
        "../scripts/tajima_d.py"



rule format_tajima_d_manhattan:
    input:
        scores=rules.tajima_d.output.scores
    output:
        formatted="results/tajima_d/tajima_d.manhattan.txt"
    log:
        "logs/format/tajima_d_manhattan.log"
    shell:
        """
        awk 'BEGIN{{OFS="\\t"}} 
             NR==1{{print "SNP", "CHR", "BP", "tajima_d"}} 
             NR>1 && $2<0 {{print "21:"$1, 21, $1, $2}}' \
        {input.scores} > {output.formatted} 2> {log}
        """


rule format_delta_tajima_d_manhattan:
    input:
        scores="results/delta_tajima_d/delta_tajima_d.scores.txt"
    output:
        formatted="results/delta_tajima_d/delta_tajima_d.manhattan.txt"
    log:
        "logs/format/delta_tajima_d_manhattan.log"
    shell:
        """
        awk 'BEGIN{{OFS="\\t"}} 
             NR==1{{print "SNP", "CHR", "BP", "delta_tajima_d", "delta_tajima_d_standardized"}} 
             NR>1 {{print "21:"$1, 21, $1, $2, $3}}' \
        {input.scores} > {output.formatted} 2> {log}
        """


rule plot_tajima_d_manhattan:
    input:
        scores="results/tajima_d/tajima_d.manhattan.txt"
    output:
        candidates="results/plots/tajima_d/tajima_d.top.candidates.txt",
        plot="results/plots/tajima_d/tajima_d.manhattan.png"
    params:
        score_column="tajima_d",
        cutoff=config["cutoff"],
        use_absolute="TRUE",
        width=640,
        height=240,
        color1="#56B4E9",
        color2="#F0E442"
    log:
        "logs/plots/tajima_d_manhattan.log"
    script:
        "../scripts/manhattan.R"


rule plot_delta_tajima_d_manhattan:
    input:
        scores="results/delta_tajima_d/delta_tajima_d.manhattan.txt"
    output:
        candidates="results/plots/delta_tajima_d/delta_tajima_d.top.candidates.txt",
        plot="results/plots/delta_tajima_d/delta_tajima_d.manhattan.png"
    params:
        score_column="delta_tajima_d_standardized",
        cutoff=config["cutoff"],
        use_absolute="TRUE",
        width=640,
        height=240,
        color1="#56B4E9",
        color2="#F0E442"
    log:
        "logs/plots/delta_tajima_d_manhattan.log"
    script:
        "../scripts/manhattan.R"
