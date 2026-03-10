using Graphs, GLMakie, GraphMakie
tree4 = [0.0 0.5 0.5 0.5;
    1.0 0.0 0.0 0.0;
    1.0 0.0 0.0 0.0;
    1.0 0.0 0.0 0.0]
intraK_4 = [
    0.000 0.409 -0.176 -0.064
    0.229 0.000 0.480 -0.404
    -0.248 0.291 0.000 -0.509
    -0.045 0.039 0.345 0.000]
tg = SimpleDiGraph(intraK_4)
fig, ax, p = graphplot(tg;
    node_size=20,
    ilabels=["1", "2", "3", "4"],
    curve_distance=0.2,  # エッジを曲線にして分離
    arrow_size=15,
    arrow_shift=:end
)
hidedecorations!(ax);
display(fig)