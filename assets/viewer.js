window.dash_clientside = window.dash_clientside || {};

window.dash_clientside.polymer_viewer = {
    installMouseControls: function (element, viewer) {
        var drag = null;
        // Capture right-button gestures before 3Dmol's default right-drag zoom.
        // Document listeners allow a drag to finish outside the viewer. Panels
        // remain mounted and this is installed only once per viewer.
        element.addEventListener("mousedown", function (event) {
            if (event.button !== 2) { return; }
            event.preventDefault();
            event.stopImmediatePropagation();
            drag = { x: event.clientX, y: event.clientY };
            element.style.cursor = "grabbing";
        }, true);
        document.addEventListener("mousemove", function (event) {
            if (!drag) { return; }
            event.preventDefault();
            event.stopImmediatePropagation();
            viewer.translateScene(event.clientX - drag.x, event.clientY - drag.y);
            drag = { x: event.clientX, y: event.clientY };
        }, true);
        function endDrag() {
            drag = null;
            element.style.cursor = "";
        }
        document.addEventListener("mouseup", function (event) {
            if (!drag) { return; }
            event.preventDefault();
            event.stopImmediatePropagation();
            endDrag();
        }, true);
        window.addEventListener("blur", endDrag);
        element.addEventListener("contextmenu", function (event) {
            event.preventDefault();
        });
    },

    // Renders a PDB block into the 3Dmol.js viewer at `elementId`. If `storeId`
    // is given, clicking an atom writes {index, elem} into that dcc.Store via
    // dash_clientside.set_props -- the bridge from JS click events back into a
    // normal server-side Dash callback, with no custom component needed.
    render: function (pdbblock, elementId, storeId, headIdx, tailIdx, clickable, uploadId, pickingRole) {
        var element = document.getElementById(elementId);
        if (!element) {
            return window.dash_clientside.no_update;
        }

        var viewer = element._3dmolViewer;
        if (!pdbblock) {
            if (viewer) {
                viewer.removeAllModels();
                viewer.render();
            }
            element._pdbblock = null;
            element._uploadId = null;
            return Date.now();
        }
        var savedView = viewer && element._pdbblock === pdbblock && element._uploadId === uploadId
            ? viewer.getView() : null;
        if (!viewer) {
            viewer = $3Dmol.createViewer(element, { backgroundColor: "white" });
            element._3dmolViewer = viewer;
            this.installMouseControls(element, viewer);
            if (typeof ResizeObserver !== "undefined") {
                element._viewerResizeObserver = new ResizeObserver(function () {
                    if (element.clientWidth && element.clientHeight) {
                        viewer.resize();
                        viewer.render();
                    }
                });
                element._viewerResizeObserver.observe(element);
            }
        }

        viewer.removeAllModels();
        // keepH: 3Dmol.js's PDB parser silently drops hydrogen atoms by default.
        viewer.addModel(pdbblock, "pdb", { keepH: true });
        viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });

        if (headIdx !== null && headIdx !== undefined) {
            viewer.setStyle({ index: headIdx }, { stick: { radius: 0.15, color: "dodgerblue" }, sphere: { scale: 0.45, color: "dodgerblue" } });
        }
        if (tailIdx !== null && tailIdx !== undefined) {
            viewer.setStyle({ index: tailIdx }, { stick: { radius: 0.15, color: "crimson" }, sphere: { scale: 0.45, color: "crimson" } });
        }

        if (clickable && storeId) {
            viewer.setClickable({}, true, function (atom, clickedViewer, event) {
                if (event && event.type.indexOf("mouse") === 0 && event.button !== 0) { return; }
                var data = { index: atom.index, elem: atom.elem };
                if (uploadId) { data.upload_id = uploadId; }
                if (pickingRole) { data.role = pickingRole; }
                window.dash_clientside.set_props(storeId, { data: data });
            });
        }

        element._pdbblock = pdbblock;
        element._uploadId = uploadId;
        if (savedView) { viewer.setView(savedView); }
        else { viewer.zoomTo(); }
        viewer.render();
        return Date.now();
    },
};
