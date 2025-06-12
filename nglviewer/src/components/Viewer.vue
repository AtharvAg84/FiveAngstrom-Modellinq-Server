<script setup>
import { ref, onMounted, onBeforeUnmount, watch, shallowRef } from "vue";
import * as NGL from "ngl";
import Setting from "./Setting.vue";
import Navigate from "./Navigate.vue";

// declarations
const stage = ref(null);
const viewer = shallowRef(null);
const structureComponent = ref(null);
const isViewerSetting = ref(false);
const representations = ref([
    { name: "backbone", visible: false },
    { name: "cartoon", visible: false },
    { name: "line", visible: true },
    { name: "ball+stick", visible: false },
    { name: "label", visible: false },
]);
k
const getThemeBackground = () => {
    const theme = document.documentElement.getAttribute("data-theme");
    return theme === "dark" ? "#1a1a1a" : "#f5f5f5";
};

const props = defineProps({
    isSidebarOpen: {
        type: Boolean,
        required: true,
    },
    sequenceData: {
        type: String,
        default: "",
    },
    processId: {
        type: String,
        default: null,
    },
    totalSamples: {
        type: Number,
        default: 0,
    },
    minimize: {
        type: Number,
        default: 0,
    },
});

// Initialize viewer
const initViewer = () => {
    if (!viewer.value && stage.value) {
        viewer.value = new NGL.Stage(stage.value);
        viewer.value.setParameters({ backgroundColor: getThemeBackground() });
    }
};

const cleanupViewer = async () => {
    if (viewer.value) {
        // Remove all components
        await viewer.value.removeAllComponents();
        structureComponent.value = null;
    }
};

const loadMolecule = async (pdbData) => {
    if (!pdbData || !viewer.value) return;

    try {
        // Clean up existing content first
        await cleanupViewer();

        const blob = new Blob([pdbData], { type: "text/plain" });
        const component = await viewer.value.loadFile(blob, {
            ext: "pdb",
        });
        structureComponent.value = component;

        representations.value.forEach((rep) => {
            if (rep.visible && structureComponent.value) {
                structureComponent.value.addRepresentation(rep.name);
            }
        });
        viewer.value.autoView();
    } catch (error) {
        console.error("Error loading molecule:", error);
    }
};

// Watch for sequence data changes
watch(
    () => props.sequenceData,
    async (newPdbData) => {
        if (newPdbData) {
            await loadMolecule(newPdbData);
        } else {
            await cleanupViewer();
        }
    },
    { immediate: true },
);

onMounted(() => {
    initViewer();

    if (props.sequenceData) {
        loadMolecule(props.sequenceData);
    }

    // Theme observer
    const observer = new MutationObserver((mutations) => {
        mutations.forEach((mutation) => {
            if (mutation.attributeName === "data-theme" && viewer.value) {
                viewer.value.setParameters({
                    backgroundColor: getThemeBackground(),
                });
            }
        });
    });

    observer.observe(document.documentElement, {
        attributes: true,
        attributeFilter: ["data-theme"],
    });
});

onBeforeUnmount(async () => {
    await cleanupViewer();
    if (viewer.value) {
        viewer.value.dispose();
        viewer.value = null;
    }
});

const updateSequence = (newSequenceData) => {
    if (newSequenceData) {
        loadMolecule(newSequenceData);
    }
};

const toggleSettings = () => {
    isViewerSetting.value = !isViewerSetting.value;
};

const handleRepresentationUpdate = (representationName, isVisible) => {
    const representation = representations.value.find(
        (rep) => rep.name === representationName,
    );

    if (representation && structureComponent.value) {
        representation.visible = isVisible;

        if (structureComponent.value.reprList) {
            const existingRepresentation =
                structureComponent.value.reprList.find(
                    (rep) => rep.name === representationName,
                );

            if (isVisible) {
                if (!existingRepresentation) {
                    structureComponent.value.addRepresentation(
                        representationName,
                    );
                }
            } else {
                if (existingRepresentation) {
                    const reps = structureComponent.value.reprList.filter(
                        (rep) => rep.name === representationName,
                    );
                    reps.forEach((rep) =>
                        structureComponent.value.removeRepresentation(rep),
                    );
                }
            }
        }
    }
};
</script>

<template>
    <div>
        <div
            id="viewport"
            style="width: 100%; height: 100vh; overflow: hidden"
            ref="stage"
        ></div>
        <button
            v-show="!isViewerSetting"
            class="viewer-settings"
            @click="toggleSettings"
        >
            <i class="fas fa-cog"></i>
        </button>
        <Setting
            v-if="isViewerSetting"
            @close-settings="toggleSettings"
            @update-visibility="handleRepresentationUpdate"
            :isSidebarOpen="isSidebarOpen"
        />
    </div>
    <Navigate
        v-if="props.processId && props.totalSamples > 0"
        :process-id="props.processId"
        :total-samples="props.totalSamples"
        @update-sequence="updateSequence"
    />
</template>

<style scoped>
.viewer-settings {
    position: fixed;
    top: 1rem;
    right: 1rem;
    font-size: 1.3em;
    cursor: pointer;
    background: none;
    border: none;
    padding: 0.5rem;
    color: var(--text-color);
    z-index: 3;
}

#viewport {
    margin: 0rem;
    z-index: 3;
    padding: 0rem;
}
</style>
