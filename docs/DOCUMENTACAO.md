# Documentação: Algoritmo de Shor com Qiskit

Atividade 04 de Computação Quântica. Esta documentação explica a teoria usada e cada função de `shor.py`. Os comentários longos que existiam no código foram movidos para cá.

## 1. Visão geral

O algoritmo de Shor fatora um inteiro composto N em tempo polinomial. Ele é **híbrido**: só a busca do período é quântica, o resto é clássico.

| Etapa | Onde roda | Função no código |
|---|---|---|
| Escolher `a` coprimo com N | Clássico | `main` (verifica `gcd(a, N)`) |
| Encontrar o período r de `a^x mod N` | **Quântico** | `build_shor_circuit` |
| Converter a medida em r | Clássico | `estima_periodo` |
| Obter os fatores com `gcd` | Clássico | `fatores_a_partir_do_periodo` |

## 2. Matemática (exemplo N = 15, a = 7)

1. Sequência `7^x mod 15`: 1, 7, 4, 13, 1, 7, 4, 13, ... O período é **r = 4**.
2. Como r é par, `7^4 ≡ 1 (mod 15)` e então `(7^2 − 1)(7^2 + 1) ≡ 0 (mod 15)`.
3. `gcd(48, 15) = 3` e `gcd(50, 15) = 5`. Logo **15 = 3 × 5**.

O método falha quando r é ímpar ou quando `a^(r/2) ≡ −1 (mod N)`. Exemplo: N = 15, a = 14 (14 ≡ −1) dá r = 2 e resultado inconclusivo. Nesse caso basta escolher outro `a`. Por isso o algoritmo é probabilístico.

## 3. Estrutura do circuito

Para N com n bits (`N.bit_length()`):

- **Registrador de contagem**: `t = 2n` qubits. Guarda o expoente x em superposição e é o que se mede no final.
- **Registrador de trabalho**: n qubits. Guarda `a^x mod N`.
- **Registrador clássico**: t bits para o resultado da medição.

Para N = 15: n = 4, t = 8, total de 12 qubits.

## 4. Funções do código

### 4.1 `build_shor_circuit(N, a, t=None)`

Monta o circuito em quatro etapas:

1. **Preparação.** Hadamard em todos os qubits de contagem (superposição uniforme de todos os x) e porta X no primeiro qubit de trabalho, que fica em |1⟩.
2. **Exponenciação modular controlada.** O qubit de contagem `k` controla a multiplicação por `a^(2^k) mod N`. O efeito total é `|x⟩|1⟩ → |x⟩|a^x mod N⟩`. Os multiplicadores são calculados com `pow(a, 2**k, N)`, e por isso, quando `a^(2^k) ≡ 1`, a porta aparece como "x1 mod N" no desenho (é a identidade).
3. **QFT inversa** no registrador de contagem, que transforma o padrão periódico em picos.
4. **Medição** dos qubits de contagem.

### 4.2 `controlled_mult_gate(mult, N, n)`

Cria uma porta controlada que faz `|y⟩ → |(mult · y) mod N⟩` no registrador de trabalho **somente se** o qubit de controle for 1.

- Como `mult` é coprimo com N, a operação é uma permutação dos estados com y < N, portanto é unitária.
- Estados com y ≥ N não são usados e ficam inalterados, para a matriz continuar unitária.
- A porta é construída como matriz de permutação de dimensão `2^(n+1)`. Convenção do Qiskit: o qubit 0 da porta é o bit menos significativo do índice, e aqui ele é o **controle**.

**Limitação (declarar no relatório):** a matriz cresce exponencialmente com n. Serve para N pequeno (15, 21, 35). No algoritmo de Shor real, essa parte é um circuito aritmético com somadores quânticos, e é a etapa mais cara.

### 4.3 `qft_dagger(t)`

QFT inversa em `t` qubits:

1. Troca (`swap`) dos qubits simétricos, para inverter a ordem dos bits.
2. Para cada qubit `j`: rotações de fase controladas `cp(−π / 2^(j−m))` vindas dos qubits anteriores `m < j`, seguidas de Hadamard.

Intuição: a QFT é o análogo quântico da transformada de Fourier. Ela detecta a "frequência" do padrão periódico do registrador de contagem. Custa cerca de t² portas, contra cerca de t·2^t operações da FFT clássica.

### 4.4 `estima_periodo(y, t, N, a)`

Converte a medida `y` em candidato a período:

1. Fase: `y / 2^t`. Com alta probabilidade ela vale aproximadamente `k / r`, com k inteiro.
2. `Fraction(fase).limit_denominator(N)` aplica frações contínuas e devolve a melhor fração com denominador menor ou igual a N. O denominador é o candidato a r.
3. Como a fração pode ter sido reduzida (ex.: 3/6 vira 1/2), testa múltiplos de r (até 5r) e fica com o primeiro que satisfaz `a^r ≡ 1 (mod N)`.
4. Se `y = 0`, devolve `None`, pois a fase 0 não carrega informação sobre r.

### 4.5 `fatores_a_partir_do_periodo(N, a, r)`

Aplica as regras de pós-processamento:

1. r ímpar: falha.
2. `a^(r/2) ≡ −1 (mod N)`: falha.
3. Senão, `p = gcd(a^(r/2) − 1, N)` e `q = gcd(a^(r/2) + 1, N)`. Se algum for 1 ou N, falha. Caso contrário, retorna `(p, q)`.

### 4.6 `caminho_evidencia(nome, extensao, carimbo)`

Monta o caminho `evidencias/AAAA-MM-DD_HH-MM-SS_nome.ext`. O carimbo de data e hora é calculado **uma vez por execução** em `main`, para que todos os arquivos de uma mesma execução compartilhem o mesmo horário. A pasta é criada se não existir.

### 4.7 `main()`

Fluxo: lê N e `a` da linha de comando (padrão 15 e 7), confere o `gcd`, calcula o período clássico só para conferência, constrói e desenha o circuito, simula com o `AerSimulator` (2048 execuções), imprime a tabela de medidas, tenta obter os fatores, salva o histograma e grava o log de texto.

## 5. Evidências (log de execuções)

Cada execução gera, em `evidencias/`:

| Arquivo | Conteúdo |
|---|---|
| `AAAA-MM-DD_HH-MM-SS_circuito_shor_N<N>_a<a>.png` | Circuito completo |
| `AAAA-MM-DD_HH-MM-SS_circuito_qft_inversa_t<t>.png` | Circuito da QFT inversa |
| `AAAA-MM-DD_HH-MM-SS_histograma_N<N>_a<a>.png` | Histograma das medidas |
| `AAAA-MM-DD_HH-MM-SS_log_N<N>_a<a>.txt` | Parâmetros, período, resultado e medidas mais frequentes |

Os arquivos nunca são sobrescritos, então a pasta funciona como histórico das execuções.

## 6. Como interpretar os resultados

- **r potência de 2 (ex.: N = 15, r = 4):** só aparecem `r` valores de y (0, 64, 128, 192), todos exatamente em múltiplos de `2^t / r`.
- **r não potência de 2 (ex.: N = 21, r = 6):** os picos se espalham, porque `2^t / r` não é inteiro. As frações contínuas são necessárias.
- **Medida que dá só divisor de r:** acontece (ex.: 1/2 quando r = 6). Por isso o código testa múltiplos e tenta outras medidas.
- **y = 0:** sempre ocorre e nunca ajuda.

## 7. Limitações

1. A exponenciação modular usa matriz de permutação (não escala, ver 4.2).
2. A execução é em simulador **sem ruído**. Em hardware real os resultados seriam mais "sujos".
3. O número de qubits (12 para N = 15) é maior que o das implementações otimizadas da literatura.

## 8. Resultados obtidos nos testes

| N | a | r | Resultado |
|---|---|---|---|
| 15 | 7 | 4 | 3 × 5 |
| 15 | 14 | 2 | inconclusivo (a^(r/2) ≡ −1) |
| 21 | 2 | 6 | 7 × 3 |
| 35 | 2 | 12 | 7 × 5 |
