# 后端部署指南 - Render

本指南将帮助你将 Lotus API 后端部署到 Render（免费）。

## 前置要求

1. GitHub 账户
2. Render 账户（免费注册：https://render.com）

## 部署步骤

### 1. 准备代码

确保以下文件已提交到 GitHub：
- `Procfile`
- `runtime.txt`
- `requirements.txt`（已包含 flask, flask-cors, gunicorn）
- `lotus/api/app.py`（已修改支持 Gunicorn）
- `.render.yaml`（可选）

### 2. 在 Render 创建 Web Service

1. 登录 Render：https://dashboard.render.com
2. 点击 **New +** → **Web Service**
3. 连接你的 GitHub 仓库
4. 选择 `Lotus` 仓库
5. 配置服务：
   - **Name**: `lotus-api`（或你喜欢的名字）
   - **Environment**: `Python 3`
   - **Region**: 选择离你最近的区域
   - **Branch**: `main`
   - **Root Directory**: （留空，使用根目录）
   - **Build Command**: 
     ```
     pip install -r requirements.txt && pip install -e .
     ```
   - **Start Command**: 
     ```
     gunicorn lotus.api.app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 120
     ```
   - **Instance Type**: `Free`（免费版）

6. 添加环境变量：
   - 点击 **Environment** 标签
   - 添加：
     - Key: `GUNICORN`
     - Value: `true`

7. 点击 **Create Web Service**

### 3. 等待部署完成

- Render 会自动开始构建和部署
- 首次部署可能需要 5-10 分钟（安装依赖）
- 构建完成后，你会看到一个 URL，例如：`https://lotus-api.onrender.com`

### 4. 测试 API

部署完成后，访问：
```
https://your-app-name.onrender.com/api/health
```

应该返回：
```json
{
  "status": "ok",
  "lotus_available": true
}
```

### 5. 更新前端配置

部署后端后，需要更新前端的 API 地址：

1. 复制你的 Render URL（例如：`https://lotus-api.onrender.com`）
2. 更新 `.github/workflows/web-demo.yml` 中的 API URL：
   ```yaml
   window.API_BASE_URL = window.API_BASE_URL || 'https://lotus-api.onrender.com';
   ```
3. 提交并推送更改
4. GitHub Actions 会自动重新部署前端

## 注意事项

### 免费版限制

- **休眠机制**：15 分钟无活动后服务会休眠
- **首次访问慢**：休眠后首次访问需要 30-60 秒唤醒
- **后续访问正常**：唤醒后响应速度正常

### 文件上传

- 上传的文件存储在临时目录（`/tmp/lotus_web_demo`）
- 服务重启后文件会丢失（这是正常的，因为使用临时存储）
- 如需持久化存储，可以考虑：
  - 升级到付费计划（有持久磁盘）
  - 使用外部存储（如 AWS S3）

### 日志查看

- 在 Render Dashboard 中点击你的服务
- 查看 **Logs** 标签可以看到实时日志

### 环境变量

如果需要配置其他环境变量：
- 在 Render Dashboard → Environment 中添加
- 例如：`UPLOAD_FOLDER`（如果需要自定义上传目录）

## 故障排除

### 构建失败

1. 检查 `requirements.txt` 是否包含所有依赖
2. 查看 Render 构建日志中的错误信息
3. 确保 Python 版本正确（3.10）

### 服务无法启动

1. 检查 `Procfile` 中的命令是否正确
2. 确保 `gunicorn` 已安装（在 requirements.txt 中）
3. 查看 Render 日志中的错误信息

### API 返回 500 错误

1. 查看 Render 日志
2. 检查是否缺少依赖包
3. 确保 `lotus` 包已正确安装（`pip install -e .`）

## 下一步

部署完后端后，继续部署前端到 GitHub Pages：
- 参考 `.github/workflows/web-demo.yml`
- 更新前端中的 API URL
- 提交代码触发自动部署

