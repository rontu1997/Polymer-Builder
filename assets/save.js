(function () {
    "use strict";
    var configurations = new Map();
    var pending = new Map();
    var active = new Set();
    var formats = {
        ".pdb": { description: "PDB structure", accept: { "chemical/x-pdb": [".pdb"] } },
        ".mol": { description: "MOL structure", accept: { "chemical/x-mdl-molfile": [".mol"] } },
        ".mol2": { description: "MOL2 structure", accept: { "chemical/x-mol2": [".mol2"] } }
    };

    function status(prefix, text, error) {
        window.dash_clientside.set_props(prefix + "-status", {
            children: text, className: error ? "error-text" : ""
        });
    }

    function isCurrent(job) {
        var current = configurations.get(job.prefix);
        return current && current.ready && current.revision === job.revision;
    }

    function finish(job) {
        clearTimeout(job.timeout);
        pending.delete(job.id);
        active.delete(job.prefix);
    }

    async function startSave(prefix) {
        var config = configurations.get(prefix);
        if (!config || !config.ready || active.has(prefix)) { return; }
        active.add(prefix);
        var job = { prefix: prefix, revision: config.revision, handle: null };
        try {
            var format = config.format;
            if (typeof window.showSaveFilePicker === "function" && window.isSecureContext) {
                status(prefix, "Choose a filename and save location…");
                // This call happens directly in the click handler, before any
                // server conversion, while browser user activation is available.
                job.handle = await window.showSaveFilePicker({
                    suggestedName: config.stem + format,
                    types: [formats[format]].concat(Object.keys(formats)
                        .filter(function (key) { return key !== format; })
                        .map(function (key) { return formats[key]; })),
                    excludeAcceptAllOption: true
                });
                var extension = /\.(pdb|mol2?)$/i.exec(job.handle.name);
                if (!extension) { throw new Error("Use a .pdb, .mol, or .mol2 filename."); }
                format = "." + extension[1].toLowerCase();
            }
            if (!isCurrent(job)) { throw new Error("The structure changed. Click Save again."); }
            job.id = window.crypto && typeof window.crypto.randomUUID === "function"
                ? window.crypto.randomUUID()
                : Date.now().toString(36) + "-" + Math.random().toString(36).slice(2);
            pending.set(job.id, job);
            job.timeout = setTimeout(function () {
                finish(job);
                status(prefix, "The export timed out. Please try Save again.", true);
            }, 120000);
            status(prefix, "Preparing " + format.slice(1).toUpperCase() + " file…");
            window.dash_clientside.set_props(prefix + "-request", {
                data: { request_id: job.id, format: format }
            });
        } catch (error) {
            finish(job);
            status(prefix, error.name === "AbortError" ? "Save cancelled." : error.message,
                error.name !== "AbortError");
        }
    }

    window.polymer_save = {
        configure: function (prefix, data, format, merged) {
            var ready = Boolean(data && (merged ? data.pdbblock : data.polymer_pdbblock));
            var previous = configurations.get(prefix);
            var filename = ready ? data.filename : "structure.pdb";
            var stem = filename.replace(/^.*[\\/]/, "").replace(/\.[^.]+$/, "");
            var config = {
                ready: ready, revision: ready ? JSON.stringify(data) : null,
                stem: stem + (merged ? "" : "_polymer"),
                format: formats[format] ? format : ".pdb"
            };
            configurations.set(prefix, config);
            if (previous && previous.revision !== config.revision) { status(prefix, ""); }
            return Date.now();
        },

        receive: async function (prefix, payload) {
            if (!payload) { return window.dash_clientside.no_update; }
            var job = pending.get(payload.request_id);
            if (!job || job.prefix !== prefix) { return window.dash_clientside.no_update; }
            clearTimeout(job.timeout);
            var writable = null;
            try {
                if (payload.error) { throw new Error(payload.error); }
                if (!isCurrent(job)) { throw new Error("The structure changed. Click Save again."); }
                if (job.handle) {
                    writable = await job.handle.createWritable();
                    if (!isCurrent(job)) { throw new Error("The structure changed. Click Save again."); }
                    await writable.write(payload.content);
                    if (!isCurrent(job)) { throw new Error("The structure changed. Click Save again."); }
                    await writable.close();
                    writable = null;
                    status(prefix, "Saved " + job.handle.name + ".");
                } else {
                    var url = URL.createObjectURL(new Blob([payload.content], {type: payload.type}));
                    var link = document.createElement("a");
                    link.href = url;
                    link.download = payload.filename;
                    document.body.appendChild(link);
                    link.click();
                    link.remove();
                    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
                    status(prefix, "Download started. To choose a folder, enable your browser’s ‘Ask where to save each file’ setting.");
                }
            } catch (error) {
                if (writable) {
                    try { await writable.abort(); } catch (_) { /* Already closed or unavailable. */ }
                }
                status(prefix, "Could not save: " + error.message, true);
            } finally {
                finish(job);
            }
            return Date.now();
        }
    };

    document.addEventListener("click", function (event) {
        var button = event.target.closest && event.target.closest("#save-button, #merge-save-button");
        if (!button || button.disabled) { return; }
        startSave(button.id.replace(/-button$/, ""));
    });
})();
