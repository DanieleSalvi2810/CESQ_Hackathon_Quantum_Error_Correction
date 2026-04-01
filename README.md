# Quantum Donuts

Beamer slides for our **QIF Challenge 1** presentation on **quantum error correction**, with a focus on **toric/surface codes**, **biased-noise decoding**, **defect-based logical qubits**, and **erasure recovery**.



<p align="center">
  <a href="slides/Quantum_Donuts.pdf"><strong>Full slides (PDF)</strong></a>
</p>

---

## 1. From stabilizers to syndrome extraction


<p align="center">
  <img src="assets/readme/tmp/slide-2.png" width="900" alt="Introductioin">
</p>

<p align="center">
  <img src="assets/readme/tmp/slide-4.png" width="900" alt="Anti-commutation">
</p>


<p align="center">
  <img src="assets/readme/tmp/slide-5.png" width="900" alt="Plaquette and cross stabilizers">
</p>

<p align="center">
  <img src="assets/readme/gifs/syndromes_pairs.gif" width="900" alt="Syndromes appear in pairs">
</p>

<p align="center">
  <img src="assets/readme/gifs/infer_error.gif" width="900" alt="Inferring the most likely error path">
</p>

---

## 2. Topology as protection

A key idea is that not all closed error paths are equivalent.  
Local, contractible loops are harmless because they are products of stabilizers.  
Global, non-contractible loops wrap around the topology and therefore implement logical operations on the encoded state. The deck uses this viewpoint to explain why the toric code stores information in global topology rather than in a single local qubit.

<p align="center">
  <img src="assets/readme/gifs/toric_loops.gif" width="900" alt="Loops on the torus">
</p>

<p align="center">
  <img src="assets/readme/tmp/slide-16.png" width="900" alt="Plaquette and cross stabilizers">
</p>

<p align="center">
  <img src="assets/readme/tmp/slide-19.png" width="900" alt="Plaquette and cross stabilizers">
</p>



<p align="center">
  <img src="assets/readme/tmp/slide-23.png" width="900" alt="Plaquette and cross stabilizers">
</p>




---

## 3. Bias-adapted decoding and asymmetric codes



<p align="center">
  <img src="assets/readme/tmp/slide-25.png" width="900" alt="Unbalanced errors in symmetric codes">
</p>

<p align="center">
  <img src="assets/readme/tmp/slide-28.png" width="900" alt="Asymmetric surface code">
</p>

<p align="center">
  <img src="assets/readme/gifs/tailored_error_chain.gif" width="900" alt="Tailored error chains">
</p>

<p align="center">
  <img src="assets/readme/gifs/biased_matching.gif" width="900" alt="Biased syndrome matching">
</p>

---

## 4. Creating additional logical qubits with defects



<p align="center">
  <img src="assets/readme/tmp/slide-38.png" width="900" alt="Creating logical qubits via defects">
</p>

<p align="center">
  <img src="assets/readme/tmp/slide-39.png" width="900" alt="Scaling larger logical qubits">
</p>

---

## 5. Erasure recovery



<p align="center">
  <img src="assets/readme/gifs/erasures_recovery.gif" width="900" alt="Recovering from erasures">
</p>

<p align="center">
  <img src="assets/readme/gifs/erasure_hardware.gif" width="900" alt="Hardware recovery of an erasure">
</p>

<p align="center">
  <img src="assets/readme/gifs/erasure_failure.gif" width="900" alt="When erasures fail">
</p>

---

## Interactive demo


<p align="center">
  <a href="https://cesq-hackathon-toric-codes-for-quantum.onrender.com/">Open the demo</a>
</p>

---

<details>
<summary>Extra static slides</summary>

<br>

<p align="center">
  <img src="assets/readme/tmp/slide-1.png" width="900" alt="Title slide">
</p>

<p align="center">
  <img src="assets/readme/tmp/slide-6.png" width="900" alt="Toric code overview">
</p>

<p align="center">
  <img src="assets/readme/tmp/slide-51.png" width="900" alt="Interactive demo QR">
</p>

</details>