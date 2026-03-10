// Web Worker for LOOM FM-index live search.
// Loads a serialized .idx file via WASM and performs exact-match search.

import init, { LoomWasmIndex } from "./pkg/loom.js";

let index = null;
let wasmReady = false;

async function ensureWasm() {
  if (!wasmReady) {
    await init();
    wasmReady = true;
  }
}

self.onmessage = async (e) => {
  const msg = e.data;

  if (msg.type === "load") {
    try {
      await ensureWasm();
      const resp = await fetch(msg.url);
      if (!resp.ok) throw new Error(`HTTP ${resp.status}: ${resp.statusText}`);
      const bytes = new Uint8Array(await resp.arrayBuffer());
      index = LoomWasmIndex.fromSerialized(bytes);
      self.postMessage({
        type: "loaded",
        files: index.fileCount(),
        corpus: index.corpusSize(),
      });
    } catch (err) {
      self.postMessage({ type: "error", error: err.message });
    }
  } else if (msg.type === "search") {
    if (!index) {
      self.postMessage({ type: "error", error: "No index loaded" });
      return;
    }
    try {
      const results = index.search(msg.pattern, msg.limit || 100, msg.ctxChars || 40);
      self.postMessage({
        type: "results",
        hits: results,
        pattern: msg.pattern,
        rid: msg.rid,
      });
    } catch (err) {
      self.postMessage({ type: "error", error: err.message });
    }
  }
};
