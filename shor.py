import sys
from fractions import Fraction
from math import gcd, pi

import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import UnitaryGate
from qiskit_aer import AerSimulator


# ---------------------------------------------------------------------------
# Etapa 2: exponenciacao modular controlada
# ---------------------------------------------------------------------------
def controlled_mult_gate(mult: int, N: int, n: int) -> UnitaryGate:
    dim_alvo = 2 ** n
    dim = 2 * dim_alvo
    U = np.zeros((dim, dim), dtype=complex)
    for c in range(2):
        for y in range(dim_alvo):
            if c == 1 and y < N:
                y_novo = (mult * y) % N
            else:
                y_novo = y
            U[c + 2 * y_novo, c + 2 * y] = 1.0
    return UnitaryGate(U, label=f"x{mult} mod {N}")


# ---------------------------------------------------------------------------
# Etapa 3: QFT inversa (implementada "na mao" para mostrar as portas)
# ---------------------------------------------------------------------------
def qft_dagger(t: int) -> QuantumCircuit:
    """QFT inversa em t qubits: SWAPs + Hadamards + fases controladas."""
    qc = QuantumCircuit(t, name="QFT†")
    for q in range(t // 2):
        qc.swap(q, t - q - 1)
    for j in range(t):
        for m in range(j):
            qc.cp(-pi / 2 ** (j - m), m, j)
        qc.h(j)
    return qc


# ---------------------------------------------------------------------------
# Circuito completo
# ---------------------------------------------------------------------------
def build_shor_circuit(N: int, a: int, t: int | None = None) -> QuantumCircuit:
    n = N.bit_length()          # qubits do registrador de trabalho
    if t is None:
        t = 2 * n               # qubits do registrador de contagem

    cont = QuantumRegister(t, "contagem")
    trab = QuantumRegister(n, "trabalho")
    cl = ClassicalRegister(t, "medida")
    qc = QuantumCircuit(cont, trab, cl)

    # Etapa 1: preparacao dos registradores
    qc.h(cont)                  # superposicao uniforme de todos os x
    qc.x(trab[0])               # registrador de trabalho = |0...01> = |1>
    qc.barrier(label="preparacao")

    # Etapa 2: exponenciacao modular controlada
    # o qubit k de contagem controla a multiplicacao por a^(2^k) mod N
    for k in range(t):
        mult = pow(a, 2 ** k, N)
        qc.append(controlled_mult_gate(mult, N, n), [cont[k], *trab])
    qc.barrier(label="a^x mod N")

    # Etapa 3: QFT inversa no registrador de contagem
    qc.append(qft_dagger(t).to_gate(label="QFT†"), cont)
    qc.barrier(label="QFT")

    # Etapa 4: medicao do registrador de contagem
    qc.measure(cont, cl)
    return qc


# ---------------------------------------------------------------------------
# Etapa 5: pos-processamento classico
# ---------------------------------------------------------------------------
def fatores_a_partir_do_periodo(N: int, a: int, r: int):
    """Retorna (p, q) se o periodo r permite achar fatores nao triviais."""
    if r % 2 != 0:
        return None
    meio = pow(a, r // 2, N)
    if meio == N - 1:           # a^(r/2) = -1 (mod N): inconclusivo
        return None
    p, q = gcd(meio - 1, N), gcd(meio + 1, N)
    if p in (1, N) or q in (1, N):
        return None
    return p, q


def estima_periodo(y: int, t: int, N: int, a: int):
    """Converte a medida y em um candidato a periodo r usando fracoes continuas."""
    fase = y / 2 ** t
    frac = Fraction(fase).limit_denominator(N)
    if y == 0:                  # fase 0 nao carrega informacao sobre r
        return fase, frac, None
    r = frac.denominator
    # se r nao for o periodo exato, tenta multiplos (r, 2r, 3r, ...)
    for k in range(1, 6):
        if r * k < N and pow(a, r * k, N) == 1:
            return fase, frac, r * k
    return fase, frac, None


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    a = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    SHOTS = 2048

    print(f"=== Algoritmo de Shor: N = {N}, a = {a} ===")
    if gcd(a, N) != 1:
        print(f"mdc({a},{N}) = {gcd(a, N)} > 1: ja e um fator, sem precisar do quantico!")
        return

    # Periodo classico (apenas para conferencia e para o relatorio)
    seq = [pow(a, x, N) for x in range(2 * N)]
    r_real = next(x for x in range(1, N + 1) if pow(a, x, N) == 1)
    print(f"Sequencia a^x mod N: {seq[:r_real * 2 + 1]} ...")
    print(f"Periodo classico (conferencia): r = {r_real}\n")

    qc = build_shor_circuit(N, a)
    t = qc.num_clbits
    print(f"Qubits: {qc.num_qubits} (contagem = {t}, trabalho = {qc.num_qubits - t})")

    # Desenho dos circuitos (evidencias para o relatorio)
    try:
        qc.draw("mpl", filename=f"circuito_shor_N{N}_a{a}.png", fold=-1)
        qft_dagger(t).draw("mpl", filename=f"circuito_qft_inversa_t{t}.png")
        print("Imagens salvas: circuito_shor / circuito_qft_inversa")
    except Exception as e:  # matplotlib/pylatexenc ausentes
        print("Nao foi possivel desenhar (pip install matplotlib pylatexenc):", e)

    # Execucao no simulador
    sim = AerSimulator()
    tqc = transpile(qc, sim)
    counts = sim.run(tqc, shots=SHOTS).result().get_counts()

    # Tabela de resultados
    print(f"\nMedidas (top 8 de {len(counts)} resultados distintos, {SHOTS} execucoes):")
    print(f"{'y (bin)':>{t + 2}} {'y':>6} {'cont.':>6} {'y/2^t':>9} {'fracao':>9} {'r':>4}")
    ordenado = sorted(counts.items(), key=lambda kv: -kv[1])
    for bits, c in ordenado[:8]:
        y = int(bits, 2)
        fase, frac, r = estima_periodo(y, t, N, a)
        print(f"{bits:>{t + 2}} {y:>6} {c:>6} {fase:>9.4f} {str(frac):>9} {str(r):>4}")

    # Estimativa do periodo e fatores
    print("\nTentando obter os fatores a partir das medidas (mais frequentes primeiro):")
    for bits, c in ordenado:
        y = int(bits, 2)
        if y == 0:
            continue            # y = 0 nao informa nada sobre r
        fase, frac, r = estima_periodo(y, t, N, a)
        if r is None:
            continue
        res = fatores_a_partir_do_periodo(N, a, r)
        print(f"  y={y}: y/2^t={fase:.4f} ~ {frac} -> r={r}; ", end="")
        if res:
            p, q = res
            print(f"mdc({a}^{r // 2}-1,{N})={p}, mdc({a}^{r // 2}+1,{N})={q}")
            print(f"\n>>> SUCESSO: {N} = {p} x {q}")
            break
        print("inconclusivo (r impar ou a^(r/2) = -1); tentando outra medida")
    else:
        print("\n>>> Nenhuma medida deu fatores: rode de novo ou escolha outro 'a'.")

    try:
        from qiskit.visualization import plot_histogram
        plot_histogram(counts, figsize=(12, 5), title=f"Medidas - N={N}, a={a}").savefig(
            f"histograma_N{N}_a{a}.png", bbox_inches="tight")
        print(f"Histograma salvo: histograma_N{N}_a{a}.png")
    except Exception as e:
        print("Histograma nao gerado:", e)


if __name__ == "__main__":
    main()