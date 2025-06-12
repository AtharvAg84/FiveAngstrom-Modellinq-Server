<script setup>
import { ref, watch } from "vue";
import Sidebar from "./components/Sidebar.vue";
import Layout from "./components/Layout.vue";
import Config from "./components/Config.vue";
import Home from "./components/Home.vue";
import Viewer from "./components/Viewer.vue";

const isSidebarOpen = ref(true);
const sequenceData = ref("");
const activeView = ref("home");
const processId = ref(null);
const totalSamples = ref(0);
const minimize = ref(0);
const key = ref(0); // Add key

const toggleSidebar = () => {
    isSidebarOpen.value = !isSidebarOpen.value;
};

// Function to update the sequence data
const updateSequenceData = (sequence) => {
    sequenceData.value = sequence;
};
// Function to change the active view
const changeView = (view) => {
    activeView.value = view;
};

const gotoHome = () => {
    activeView.value = "home";
};

const updateProcessInfo = (processInfo) => {
    processId.value = processInfo.processId;
    totalSamples.value = processInfo.totalSamples;
    minimize.value = processInfo.minimize;
    key.value += 1;
};
</script>

<template>
    <Layout :isSidebarOpen="isSidebarOpen" @toggle-sidebar="toggleSidebar">
        <Sidebar
            :is-sidebar-open="isSidebarOpen"
            @toggle-sidebar="toggleSidebar"
            @goto-home="gotoHome"
        >
            <Config
                @sequence-data="updateSequenceData"
                @process-info="updateProcessInfo"
            />
        </Sidebar>
        <Home
            v-if="activeView === 'home'"
            :is-sidebar-open="isSidebarOpen"
            @start-viewer="changeView('viewer')"
        />
        <Viewer
            v-if="activeView === 'viewer'"
            :key="key"
            :isSidebarOpen="isSidebarOpen"
            :sequenceData="sequenceData"
            :processId="processId"
            :totalSamples="totalSamples"
            :minimize="minimize"
        />
    </Layout>
</template>
