import sys
from datetime import datetime
from fractions import Fraction
from math import gcd, pi
from pathlib import Path

import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit_aer import AerSimulator

PASTA_EVIDENCIAS = Path(__file__).resolve().parent / "evidencias"
SHOTS = 2048


def caminho_evidencia(nome: str, extensao: str, carimbo: str) -> Path:
    """Caminho em evidencias/ no formato AAAA-MM-DD_HH-MM-SS_nome.ext."""
    PASTA_EVIDENCIAS.mkdir(exist_ok=True)
    return PASTA_EVIDENCIAS / f"{carimbo}_{nome}.{extensao}"


def controlled_mult_gate(mult: int, N: int, n: int) -> UnitaryGate:
    """Porta controlada |y> -> |(mult*y) mod N> (ver docs, secao 4.2)."""
    dim_alvo = 2 ** n
    U = np.zeros((2 * dim_alvo, 2 * dim_alvo), dtype=complex)
    for c in range(2):
        for y in range(dim_alvo):
            y_novo = (mult * y) % N if (c == 1 and y < N) else y
            U[c + 2 * y_novo, c + 2 * y] = 1.0
    return UnitaryGate(U, label=f"x{mult} mod {N}")


def qft_dagger(t: int) -> QuantumCircuit:
    """QFT inversa em t qubits (ver docs, secao 4.3)."""
    qc = QuantumCircuit(t, name="QFT†")
    for q in range(t // 2):
        qc.swap(q, t - q - 1)
    for j in range(t):
        for m in range(j):
            qc.cp(-pi / 2 ** (j - m), m, j)
        qc.h(j)
    return qc


def build_shor_circuit(N: int, a: int, t: int | None = None) -> QuantumCircuit:
    """Monta o circuito completo (ver docs, secao 4.1)."""
    n = N.bit_length()
    t = 2 * n if t is None else t

    cont = QuantumRegister(t, "contagem")
    trab = QuantumRegister(n, "trabalho")
    cl = ClassicalRegister(t, "medida")
    qc = QuantumCircuit(cont, trab, cl)

    qc.h(cont)
    qc.x(trab[0])
    qc.barrier(label="preparacao")

    for k in range(t):
        qc.append(controlled_mult_gate(pow(a, 2 ** k, N), N, n), [cont[k], *trab])
    qc.barrier(label="a^x mod N")

    qc.append(qft_dagger(t).to_gate(label="QFT†"), cont)
    qc.barrier(label="QFT")

    qc.measure(cont, cl)
    return qc


def fatores_a_partir_do_periodo(N: int, a: int, r: int):
    """Retorna (p, q) se r gerar fatores nao triviais; senao None."""
    if r % 2 != 0:
        return None
    meio = pow(a, r // 2, N)
    if meio == N - 1:
        return None
    p, q = gcd(meio - 1, N), gcd(meio + 1, N)
    if p in (1, N) or q in (1, N):
        return None
    return p, q


def estima_periodo(y: int, t: int, N: int, a: int):
    """Medida y -> (fase, fracao, r) por fracoes continuas (docs, secao 4.4)."""
    fase = y / 2 ** t
    frac = Fraction(fase).limit_denominator(N)
    if y == 0:
        return fase, frac, None
    r = frac.denominator
    for k in range(1, 6):
        if r * k < N and pow(a, r * k, N) == 1:
            return fase, frac, r * k
    return fase, frac, None


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    a = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    carimbo = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    print(f"=== Algoritmo de Shor: N = {N}, a = {a} ({carimbo}) ===")
    if gcd(a, N) != 1:
        print(f"mdc({a},{N}) = {gcd(a, N)} > 1: ja e um fator, sem precisar do quantico!")
        return

    r_real = next(x for x in range(1, N + 1) if pow(a, x, N) == 1)
    seq = [pow(a, x, N) for x in range(2 * r_real + 1)]
    print(f"Sequencia a^x mod N: {seq} ...")
    print(f"Periodo classico (conferencia): r = {r_real}\n")

    qc = build_shor_circuit(N, a)
    t = qc.num_clbits
    print(f"Qubits: {qc.num_qubits} (contagem = {t}, trabalho = {qc.num_qubits - t})")

    sufixo = f"N{N}_a{a}"
    try:
        f1 = caminho_evidencia(f"circuito_shor_{sufixo}", "png", carimbo)
        f2 = caminho_evidencia(f"circuito_qft_inversa_t{t}", "png", carimbo)
        qc.draw("mpl", filename=str(f1), fold=-1)
        qft_dagger(t).draw("mpl", filename=str(f2))
        print(f"Evidencias salvas: {f1.name}, {f2.name}")
    except Exception as e:
        print("Nao foi possivel desenhar (pip install matplotlib pylatexenc):", e)

    sim = AerSimulator()
    counts = sim.run(transpile(qc, sim), shots=SHOTS).result().get_counts()

    print(f"\nMedidas (top 8 de {len(counts)} resultados distintos, {SHOTS} execucoes):")
    print(f"{'y (bin)':>{t + 2}} {'y':>6} {'cont.':>6} {'y/2^t':>9} {'fracao':>9} {'r':>4}")
    ordenado = sorted(counts.items(), key=lambda kv: -kv[1])
    for bits, c in ordenado[:8]:
        y = int(bits, 2)
        fase, frac, r = estima_periodo(y, t, N, a)
        print(f"{bits:>{t + 2}} {y:>6} {c:>6} {fase:>9.4f} {str(frac):>9} {str(r):>4}")

    print("\nTentando obter os fatores (medidas mais frequentes primeiro):")
    resultado = None
    for bits, _ in ordenado:
        y = int(bits, 2)
        fase, frac, r = estima_periodo(y, t, N, a)
        if r is None:
            continue
        res = fatores_a_partir_do_periodo(N, a, r)
        print(f"  y={y}: y/2^t={fase:.4f} ~ {frac} -> r={r}; ", end="")
        if res:
            p, q = res
            print(f"mdc({a}^{r // 2}-1,{N})={p}, mdc({a}^{r // 2}+1,{N})={q}")
            print(f"\n>>> SUCESSO: {N} = {p} x {q}")
            resultado = res
            break
        print("inconclusivo (r impar ou a^(r/2) = -1); tentando outra medida")
    if resultado is None:
        print("\n>>> Nenhuma medida deu fatores: rode de novo ou escolha outro 'a'.")

    try:
        from qiskit.visualization import plot_histogram
        f3 = caminho_evidencia(f"histograma_{sufixo}", "png", carimbo)
        plot_histogram(counts, figsize=(12, 5), title=f"Medidas - {sufixo}").savefig(
            f3, bbox_inches="tight")
        print(f"Evidencia salva: {f3.name}")
    except Exception as e:
        print("Histograma nao gerado:", e)

    # Log textual da execucao, junto com as imagens
    log = caminho_evidencia(f"log_{sufixo}", "txt", carimbo)
    linhas = [
        f"Execucao: {carimbo}",
        f"N = {N}, a = {a}, shots = {SHOTS}, qubits = {qc.num_qubits}",
        f"Periodo classico (conferencia): r = {r_real}",
        f"Resultado: {N} = {resultado[0]} x {resultado[1]}" if resultado
        else "Resultado: inconclusivo",
        "Medidas (y: contagem): "
        + ", ".join(f"{int(b, 2)}: {c}" for b, c in ordenado[:10]),
    ]
    log.write_text("\n".join(linhas) + "\n", encoding="utf-8")
    print(f"Log salvo: {log.name}")


if __name__ == "__main__":
    main()