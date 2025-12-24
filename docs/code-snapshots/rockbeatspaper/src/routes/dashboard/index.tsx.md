---
source: /Users/byron/projects/rockbeatspaper/src/routes/dashboard/index.tsx
relative: rockbeatspaper/src/routes/dashboard/index.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { createFileRoute, Link } from "@tanstack/react-router";
import { useQuery } from "@tanstack/react-query";
import {
  Zap,
  TrendingDown,
  Thermometer,
  Activity,
  Plus,
  ArrowRight,
} from "lucide-react";
import { ProtectedLayout } from "~/components/ProtectedLayout";
import { Card } from "~/components/Card";
import { Button } from "~/components/Button";
import { useTRPC } from "~/trpc/react";
import { useAuthStore } from "~/stores/auth";

export const Route = createFileRoute("/dashboard/")({
  component: Dashboard,
});

function Dashboard() {
  const trpc = useTRPC();
  const { authToken } = useAuthStore();

  const devicesQuery = useQuery(
    trpc.getDevices.queryOptions(
      { authToken: authToken! },
      { enabled: !!authToken }
    )
  );

  const devices = devicesQuery.data?.devices || [];

  // Mock data for demonstration
  const stats = [
    {
      label: "Energy Saved",
      value: "24.5%",
      icon: TrendingDown,
      color: "from-green-600 to-teal-600",
    },
    {
      label: "Active Devices",
      value: devices.length.toString(),
      icon: Zap,
      color: "from-blue-600 to-purple-600",
    },
    {
      label: "Avg Temperature",
      value: "22°C",
      icon: Thermometer,
      color: "from-orange-600 to-red-600",
    },
    {
      label: "System Health",
      value: "98%",
      icon: Activity,
      color: "from-purple-600 to-pink-600",
    },
  ];

  return (
    <ProtectedLayout>
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {/* Header */}
        <div className="mb-8">
          <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-2">
            Dashboard
          </h1>
          <p className="text-gray-600 dark:text-gray-400">
            Monitor your energy consumption and optimize efficiency
          </p>
        </div>

        {/* Stats Grid */}
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 mb-8">
          {stats.map((stat) => {
            const Icon = stat.icon;
            return (
              <Card key={stat.label}>
                <div className="flex items-center justify-between">
                  <div>
                    <p className="text-sm text-gray-600 dark:text-gray-400 mb-1">
                      {stat.label}
                    </p>
                    <p className="text-3xl font-bold text-gray-900 dark:text-white">
                      {stat.value}
                    </p>
                  </div>
                  <div
                    className={`w-12 h-12 bg-gradient-to-r ${stat.color} rounded-lg flex items-center justify-center`}
                  >
                    <Icon className="w-6 h-6 text-white" />
                  </div>
                </div>
              </Card>
            );
          })}
        </div>

        {/* Quick Actions */}
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-8">
          <Card>
            <h2 className="text-xl font-bold text-gray-900 dark:text-white mb-4">
              AI Analysis
            </h2>
            <p className="text-gray-600 dark:text-gray-400 mb-6">
              Get intelligent recommendations to optimize your thermal systems
              and reduce energy costs.
            </p>
            <Link to="/analysis">
              <Button className="w-full sm:w-auto">
                Start Analysis
                <ArrowRight className="w-4 h-4 ml-2" />
              </Button>
            </Link>
          </Card>

          <Card>
            <h2 className="text-xl font-bold text-gray-900 dark:text-white mb-4">
              Device Management
            </h2>
            <p className="text-gray-600 dark:text-gray-400 mb-6">
              Add and manage your thermal devices to track performance and
              energy usage.
            </p>
            <Link to="/devices">
              <Button variant="outline" className="w-full sm:w-auto">
                Manage Devices
                <Plus className="w-4 h-4 ml-2" />
              </Button>
            </Link>
          </Card>
        </div>

        {/* Recent Devices */}
        <Card>
          <div className="flex items-center justify-between mb-6">
            <h2 className="text-xl font-bold text-gray-900 dark:text-white">
              Your Devices
            </h2>
            <Link to="/devices">
              <Button variant="outline" size="sm">
                View All
              </Button>
            </Link>
          </div>

          {devicesQuery.isLoading ? (
            <p className="text-gray-600 dark:text-gray-400">Loading devices...</p>
          ) : devices.length === 0 ? (
            <div className="text-center py-12">
              <Zap className="w-12 h-12 text-gray-400 mx-auto mb-4" />
              <p className="text-gray-600 dark:text-gray-400 mb-4">
                No devices yet. Add your first device to get started!
              </p>
              <Link to="/devices">
                <Button>Add Device</Button>
              </Link>
            </div>
          ) : (
            <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
              {devices.slice(0, 6).map((device) => (
                <div
                  key={device.id}
                  className="border border-gray-200 dark:border-gray-700 rounded-lg p-4 hover:shadow-md transition-shadow"
                >
                  <div className="flex items-center justify-between mb-2">
                    <h3 className="font-semibold text-gray-900 dark:text-white">
                      {device.name}
                    </h3>
                    <span
                      className={`px-2 py-1 text-xs rounded-full ${
                        device.status === "active"
                          ? "bg-green-100 text-green-800"
                          : "bg-gray-100 text-gray-800"
                      }`}
                    >
                      {device.status}
                    </span>
                  </div>
                  <p className="text-sm text-gray-600 dark:text-gray-400">
                    {device.type} • {device.location}
                  </p>
                </div>
              ))}
            </div>
          )}
        </Card>
      </div>
    </ProtectedLayout>
  );
}

```
